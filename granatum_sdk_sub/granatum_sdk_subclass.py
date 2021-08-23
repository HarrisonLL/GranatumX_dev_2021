from granatum_sdk import Granatum
import base64
from io import StringIO, BytesIO
import gzip
import json
import numpy as np
import time

class granatum_extended(Granatum):

    """
    chunk data model:
    {   "origin data size": [#, #]
        "current chunk": ["last_gbox_name", "col/row", "dense"],
        "suggested chunk": {
                "deepimpute": ["col", ##],
                "log-transform":["row", ##],
                ..
        },
        "chunk1": base64string,
        "chunk2": base64string,
        ...
    }
    """

    def __init__(self, gbox_name):
        super().__init__()
        self.gbox_name = gbox_name
        self.new_file = dict()


    def check_chunk(self, assay):
        return "current chunk" in assay


    def adjust_transform_helper2(self, chunk, count, flag, l=None, r=None):
        if flag == "col":
             new_assay = {
                     "matrix":np.array(chunk["matrix"])[:,l:r].tolist(),
                     "geneIds":chunk["geneIds"],
                     "sampleIds":chunk["sampleIds"][l:r]
                           }
             if "chunk" + str(count) in self.new_file:
                 self.new_file["chunk"+str(count)] += [self.compress_chunk(new_assay)]
             else:
                 self.new_file["chunk"+str(count)] = [self.compress_chunk(new_assay)]
        else:
             new_assay = {
                    "matrix":chunk["matrix"][l:r],
                    "geneIds":chunk["geneIds"][l:r],
                    "sampleIds":chunk["sampleIds"]
                         }
             if "chunk" + str(count) in self.new_file:
                 self.new_file["chunk"+str(count)] += [self.compress_chunk(new_assay)]
             else:
                 self.new_file["chunk"+str(count)] = [self.compress_chunk(new_assay)]


    def adjust_transform_helper1(self,assay,new_size, org_size, flag):
        if new_size < org_size:
            left = 0
            count = 1
            begin = time.time()
            for i in range(len(assay)-3):
                chunk = self.decompress_chunk(assay["chunk"+str(i+1)])
                if flag == "col":
                    current_size = len(chunk["sampleIds"])
                else:
                    current_size = len(chunk["geneIds"])
    
                if left != 0:
                     self.adjust_transform_helper2(chunk, count, flag, r=left)
                     count += 1
                     end = time.time()
                     print("finish one new chunk, it took %.3f"%(end-begin), flush=True)
                     begin = time.time()

                while True:
                    tmp = new_size * count % current_size
                    if  tmp >= new_size:
                        self.adjust_transform_helper2(chunk, count, flag, l=left, r=tmp)
                        count += 1
                        left = tmp
                        end = time.time()
                        print("finish one new chunk, it took %.3f"%(end-begin), flush=True)
                        begin = time.time()
                    elif tmp == 0:
                        self.adjust_transform_helper2(chunk, count, flag, l=left, r=current_size)
                        count += 1
                        left = tmp
                        end = time.time()
                        print("finish one new chunk, it took %.3f"%(end-begin), flush=True)
                        begin = time.time()
                        break
                    else:
                        self.adjust_transform_helper2(chunk, count, flag, l=left, r=current_size)
                        left = tmp
                        break
        else:
            count = 1
            left = new_size
            begin = time.time()
            for i in range(len(assay)-3):
                if i < len(assay) - 4:
                    current_size = org_size
                else:
                    if flag == "col":
                        current_size = len(self.decompress_chunk(assay["chunk"+str(i+1)])["sampleIds"])
                    else:
                        current_size = len(self.decompress_chunk(assay["chunk"+str(i+1)])["geneIds"])

                if  left >= current_size:
                    if "chunk"+str(count) in self.new_file:
                        self.new_file["chunk"+str(count)] += [assay["chunk"+str(i+1)]]
                    else:
                        self.new_file["chunk"+str(count)] = [assay["chunk"+str(i+1)]]
                    left -= current_size
                    if left == 0:
                        left = new_size
                        count += 1
                else:
                    chunk = self.decompress_chunk(assay["chunk"+str(i+1)])
                    self.adjust_transform_helper2(chunk, count, flag, r=left)
                    self.adjust_transform_helper2(chunk, count+1, flag, l=left)
                    left = new_size - (current_size - left)
                    count += 1
                    end = time.time()
                    print("finish one new chunk, it took %.3f"%(end-begin), flush=True)
                    begin = time.time()    


    def adjust_transform(self, assay):
        # user rerun the gbox
        if assay["current chunk"][0] == self.gbox_name:
            return
        # four situations
        if assay["current chunk"][1] == assay["suggested chunk"].get(self.gbox_name)[0] == "col":
            new_col_size = assay["suggested chunk"].get(self.gbox_name)[1]
            #uncomment the following to test
            #new_col_size = 300
            org_col_size = assay["suggested chunk"].get(assay["current chunk"][0])[1]
            self.adjust_transform_helper1(assay, new_col_size, org_col_size, "col")

        elif assay["current chunk"][1] == assay["suggested chunk"].get(self.gbox_name)[0] == "row":
            new_row_size = assay["suggested chunk"].get(self.gbox_name)[1]
            org_row_size = assay["suggested chunk"].get(assay["current chunk"][0])[1]
            self.adjust_transform_helper1(assay, new_row_size, org_row_size, "row")

        elif assay["current chunk"][1] == "col" and assay["suggested chunk"].get(self.gbox_name)[0] == "row":
            org_num_row = assay["origin data size"][0]
            sug_num_row = assay["suggested chunk"].get(self.gbox_name)[1]
            for i in range(len(assay)-3):
               chunk = self.decompress_chunk(assay["chunk"+str(i+1)])
               count = 0
               for j in range(0, org_num_row, sug_num_row):
                   count += 1
                   self.adjust_transform_helper2(chunk, count, "row", l=j, r=j+sug_num_row)

        elif assay["current chunk"][1] == "row" and assay["suggested chunk"].get(self.gbox_name)[0] == "col":
            org_num_col = assay["origin data size"][0]
            sug_num_col = assay["suggested chunk"].get(self.gbox_name)[1]
            for i in range(len(assay)-3):
               chunk = self.decompress_chunk(assay["chunk"+str(i+1)])
               count = 0
               for j in range(0, org_num_col, sug_num_col):
                   count += 1
                   self.adjust_transform_helper2(chunk, count, "col", l=j, r=j+sug_num_col)
   
    def combine_new_chunk(self, new_chunk):
        # new chunk is a list of base64string
        start_part = self.decompress_chunk(new_chunk[0])
        matrix = np.array(start_part["matrix"])
        geneIds = start_part["geneIds"]
        sampleIds = start_part["sampleIds"]
        for i in range(1, len(new_chunk)):
            part = self.decompress_chunk(new_chunk[i])
            sampleIds += part["sampleIds"]
            matrix = np.c_[matrix, np.array(part["matrix"])]
        combined_chunk = {
                "matrix":matrix.tolist(),
                "geneIds":geneIds,
                "sampleIds":sampleIds
        }
        return combined_chunk


    def decompress_chunk(self, pram1):
        bio = BytesIO()
        stream = BytesIO(base64.b64decode(pram1.encode('utf-8')))
        decompressor = gzip.GzipFile(fileobj=stream, mode='r')
        while True:  # until EOF
            chunk = decompressor.read(8192)
            if not chunk:
                decompressor.close()
                bio.seek(0)
                output = bio.read().decode("utf-8")
                break
            bio.write(chunk)
        return json.loads(output)


    def compress_chunk(self, pram1):
        tmp = json.dumps(pram1)
        bio = BytesIO()
        bio.write(tmp.encode('utf-8'))
        bio.seek(0)
        stream = BytesIO()
        compressor = gzip.GzipFile(fileobj=stream, mode = 'w')
        while True:
            chunk = bio.read(8192)
            if not chunk:
                compressor.close()
                break
            compressor.write(chunk)
        encoded = base64.b64encode(stream.getvalue())
        return(encoded.decode('utf-8'))
