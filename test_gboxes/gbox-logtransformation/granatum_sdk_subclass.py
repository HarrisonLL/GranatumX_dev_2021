from granatum_sdk import Granatum
import base64
from io import StringIO, BytesIO
import gzip
import json

class granatum_extended(Granatum):

    """
    chunk data model:
    {   "origin data size": [#, #]
        "current chunk": ["last_gbox_name", "col/row"],
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


    def check_chunk(self, assay):
        return "current chunk" in assay


    def adjust_transform(self, assay):
        # user rerun the gbox
        if assay["current chunk"][0] == self.gbox_name:
            return
        # four situations
        if assay["current chunk"][-1] == assay["suggested chunk"].get(self.gbox_name)[0] == "col":
            pass
        elif assay["current chunk"][-1] == assay["suggested chunk"].get(self.gbox_name)[0] == "row":
            pass
        elif assay["current chunk"][-1] == "col" and assay["suggested chunk"].get(self.gbox_name)[0] == "row":
            new_file = dict()
            for i in range(len(assay)-3):
               chunk = self.decompress_chunk(assay["chunk"+str(i)])
               org_num_row = assay["origin data size"][0]
               sug_num_row = assay["suggested chunk"].get(self.gbox_name)[1]
               count = 0
               for j in range(0, org_num_row, sug_num_row):
                   new_assay = {
                                "matrix":chunk["matrix"][j:j+sug_num_row],
                                "geneIds":chunk["geneIds"][j:j+sug_num_row],
                                "sampleIds":chunk["sampleIds"]
                           }
                   count += 1
                   if "chunk"+str(count) in new_file:
                       new_file["chunk"+str(count)] += [self.compress_chunk(new_assay)]
                   else:
                       new_file["chunk"+str(count)] = [self.compress_chunk(new_assay)]
            return new_file

        elif json_dict["current chunk"][-1] == "row" and json_dict["suggested chunk"].get(self.gbox_name)[0] == "col":
            new_file = dict()
            for i in range(len(assay)-3):
               chunk = self.decompress_chunk(assay["chunk"+str(i)])
               org_num_col = assay["origin data size"][0]
               sug_num_col = assay["suggested chunk"].get(self.gbox_name)[1]
               count = 0
               for j in range(0, org_num_col, sug_num_col):
                   new_assay = {
                                "matrix":np.array(chunk["matrix"][:, j:j+sug_num_col]).tolist(),
                                "geneIds":chunk["geneIds"],
                                "sampleIds":chunk["sampleIds"][j:j+sug_num_col]
                           }
                   count += 1
                   if "chunk"+str(count) in new_file:
                       new_file["chunk"+str(count)] += [self.compress_chunk(new_assay)]
                   else:
                       new_file["chunk"+str(count)] = [self.compress_chunk(new_assay)]
            return new_file


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
