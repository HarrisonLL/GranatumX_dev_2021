from granatum_sdk import Granatum
import scipy
import json
import numpy as np
import time

class granatum_extended2(Granatum):

    """
    chunk data model:
    {   "origin data size": [#, #]
        "current chunk": ["last_gbox_name", "col/row", "sparse"],
        "suggested chunk": {
                "deepimpute": ["col", ##],
                "log-transform":["row", ##],
                ..
        },
        "chunk1": {
        "sampleIds":list,
        "geneIds":list,
        "data":list,
        "indices":list,
        "indptr": list
        },
        "chunk2": {
        "sampleIds":list,
        "geneIds":list,
        "data":list,
        "indices":list,
        "indptr": list
        },
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
            row_size = len(chunk.get("geneIds"))
            col_size = len(chunk.get("sampleIds"))
            sparse = scipy.sparse.csc_matrix((chunk.get("data"), chunk.get("indices"), chunk.get("indptr")), shape=(row_size,col_size))
            tmp = sparse[:, l:r]
            new_assay = {
                "data":tmp.data,
                "indices":tmp.indices,
                "indptr":tmp.indptr,
                "geneIds":chunk["geneIds"],
                "sampleIds":chunk["sampleIds"][l:r]
                }

        else:
            row_size = len(chunk.get("geneIds"))
            col_size = len(chunk.get("sampleIds"))
            sparse = scipy.sparse.csr_matrix((chunk.get("data"), chunk.get("indices"), chunk.get("indptr")), shape=(row_size,col_size))
            tmp = sparse[l:r, :]
            new_assay = {
                "data":tmp.data,
                "indices":tmp.indices,
                "indptr":tmp.indptr,
                "geneIds":chunk["geneIds"][l:r],
                "sampleIds":chunk["sampleIds"]
                }
        
        if "chunk" + str(count) in self.new_file:
            self.new_file["chunk"+str(count)] += [new_assay]
        else:
            self.new_file["chunk"+str(count)] = [new_assay]


    def adjust_transform_helper1(self,assay,new_size, org_size, flag):
        if new_size < org_size:
            left = 0
            count = 1
            begin = time.time()
            for i in range(len(assay)-3):
                chunk = assay["chunk"+str(i+1)]
                if flag == "col":
                    current_size = len(chunk["sampleIds"])
                else:
                    current_size = len(chunk["geneIds"])
    
                if left != 0:
                     self.adjust_transform_helper2(chunk, count, flag, r=left)
                     count += 1

                while True:
                    tmp = new_size * count % current_size
                    if  tmp >= new_size:
                        self.adjust_transform_helper2(chunk, count, flag, l=left, r=tmp)
                        count += 1
                        left = tmp
                    elif tmp == 0:
                        self.adjust_transform_helper2(chunk, count, flag, l=left, r=current_size)
                        count += 1
                        left = tmp
                        break
                    else:
                        self.adjust_transform_helper2(chunk, count, flag, l=left, r=current_size)
                        left = tmp
                        break
        else:
            count = 1
            left = new_size
            for i in range(len(assay)-3):
                if i < len(assay) - 4:
                    current_size = org_size
                else:
                    if flag == "col":
                        current_size = len((assay["chunk"+str(i+1)])["sampleIds"])
                    else:
                        current_size = len((assay["chunk"+str(i+1)])["geneIds"])

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
                    chunk = assay["chunk"+str(i+1)]
                    self.adjust_transform_helper2(chunk, count, flag, r=left)
                    self.adjust_transform_helper2(chunk, count+1, flag, l=left)
                    left = new_size - (current_size - left)
                    count += 1  


    def adjust_transform(self, assay):
        # DELETE this when test is complete
        assay["suggested chunk"].update({"gbox-scanpygenefilering":["row", 500]})
        # user rerun the gbox
        if assay["current chunk"][0] == self.gbox_name:
            return
        # four situations
        if assay["current chunk"][1] == assay["suggested chunk"].get(self.gbox_name)[0] == "col":
            new_col_size = assay["suggested chunk"].get(self.gbox_name)[1]
            #uncomment the following to test
            #new_col_size = 300
            org_col_size = assay["suggested chunk"].get(assay["current chunk"][0])[1]
            if type(org_col_size) == str: org_col_size = int(org_col_size)
            if type(new_col_size) == str: new_col_size = int(new_col_size)
            self.adjust_transform_helper1(assay, new_col_size, org_col_size, "col")

        elif assay["current chunk"][1] == assay["suggested chunk"].get(self.gbox_name)[0] == "row":
            new_row_size = assay["suggested chunk"].get(self.gbox_name)[1]
            org_row_size = assay["suggested chunk"].get(assay["current chunk"][0])[1]
            if type(org_row_size) == str: org_row_size = int(org_row_size)
            if type(new_row_size) == str: new_row_size = int(new_row_size)
            self.adjust_transform_helper1(assay, new_row_size, org_row_size, "row")

        elif assay["current chunk"][1] == "col" and assay["suggested chunk"].get(self.gbox_name)[0] == "row":

            row = assay["suggested chunk"].get(self.gbox_name)[1]
            for i in range(len(assay)-3):
                chunk = assay["chunk"+str(i+1)]
                row_size = len(chunk.get("geneIds"))
                col_size = len(chunk.get("sampleIds"))
                sparse = scipy.sparse.csc_matrix((chunk.get("data"), chunk.get("indices"), chunk.get("indptr")), shape=(row_size,col_size))
                sparse_csr = sparse.tocsr()
                count = 1
                for i in range(0,row_size,row):
                    if i+row >= row_size:
                        tmp = sparse_csr[i:row_size, :]
                    else:
                        tmp = sparse_csr[i:i+row, :]
                    new_chunk = {
                        "sampleIds":chunk.get("sampleIds"),
                        "geneIds":chunk.get("geneIds")[i:i+row],
                        "data":tmp.data,
                        "indices":tmp.indices,
                        "indptr":tmp.indptr
                    }
                    if "chunk"+str(count) in self.new_file:
                        self.new_file["chunk"+str(count)] += [new_chunk]
                    else:
                        self.new_file["chunk"+str(count)] = [new_chunk]
                    count += 1

        elif assay["current chunk"][1] == "row" and assay["suggested chunk"].get(self.gbox_name)[0] == "col":
            
            column = assay["suggested chunk"].get(self.gbox_name)[1]
            for i in range(len(assay)-3):
                chunk = assay["chunk"+str(i+1)]
                row_size = len(chunk.get("geneIds"))
                col_size = len(chunk.get("sampleIds"))
                sparse = scipy.sparse.csr_matrix((chunk.get("data"), chunk.get("indices"), chunk.get("indptr")), shape=(row_size,col_size))
                sparse_csr = sparse.tocsc()
                count = 1
                for i in range(0,col_size,column):
                    if i+column >= col_size:
                        tmp = sparse_csr[:, i:col_size]
                    else:
                        tmp = sparse_csr[:, i:i+column]
                    new_chunk = {
                        "sampleIds":chunk.get("sampleIds")[i:i+column],
                        "geneIds":chunk.get("geneIds"),
                        "data":tmp.data,
                        "indices":tmp.indices,
                        "indptr":tmp.indptr
                    }
                    if "chunk"+str(count) in self.new_file:
                        self.new_file["chunk"+str(count)] += [new_chunk]
                    else:
                        self.new_file["chunk"+str(count)] = [new_chunk]
                    count += 1

   # Note the flag is the original chunk kind
    def combine_new_chunk(self, new_chunk, flag):
        sampleIds = []
        geneIds = []
        lst = []
        for piece in new_chunk:
            sampleIds += piece["sampleIds"]
            geneIds = piece["geneIds"]
            if flag == "col":
                tmp = scipy.sparse.csc_matrix((piece.get("data"), piece.get("indices"), piece.get("indptr")), shape=(len(piece["geneIds"]), len(piece["sampleIds"])))
            else:
                tmp = scipy.sparse.csr_matrix((piece.get("data"), piece.get("indices"), piece.get("indptr")), shape=(len(piece["geneIds"]), len(piece["sampleIds"])))
            lst.append(tmp)
        if flag == "col":
            combined = scipy.sparse.hstack(tuple(lst))
        else:
            combined = scipy.sparse.vstack(tuple(lst))
        
        new_chunk = {
            "sampleIds":sampleIds,
            "geneIds":geneIds,
            "data":combined.data,
            "indices":combined.indices,
            "indptr":combined.indptr
    
                }
        return new_chunk

