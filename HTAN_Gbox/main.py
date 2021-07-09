import granatum_sdk
import os
import sys
from tqdm import tqdm
import shutil
import copy
import synapseclient
import pandas as pd
import numpy as np
import scipy.io
import scipy.sparse
import gzip

def download_data(SIDs, NUM):
    # destroy if exits, and then create ./tmp_datasets directory
    dirpath = './tmp_datasets'
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
    os.mkdir(dirpath)

    syn = synapseclient.Synapse()
    syn.login("HarrisL", "Hlld@0217")

    SID_list = SIDs.replace(' ', '').split(";")

    if len(SID_list) != NUM:
        print("Error input, the number of IDs doesn't consistent with ID input.", file=sys.stderr)

    # add try except here for invalid ID

    Paths = []

    try:
        for i in range(NUM):
            entity = syn.get(SID_list[i], downloadLocation=dirpath)
            Paths.append(entity.path)
            print("Download barcode file path is:" + entity.path, flush=True)
    except ValueError:
        print("Invalid Synapse ID, please double check and copy the ID directly from the website.")

    return Paths
    
def read_files(file_name):
    file_form = file_name.split(".")[-1]
    label = ''.join(file_name.split(".")[:-1])
    print(file_form, flush = True)
    print(label, flush = True)
    if file_form == "mtx":
        Matrix = (scipy.io.mmread(file_name))
        gene_nums = Matrix.shape[0]
        cell_nums = 10000
        data = np.zeros((gene_nums, cell_nums))
        cols = Matrix.col.tolist()
        rows = Matrix.row.tolist()
        values = Matrix.data.tolist()
        for i in range(len(cols)):
            if cols[i] < cell_nums:
                data[rows[i]][cols[i]] = values[i]
                #if i % 1000 == 0:
                #    print (data[rows[i]][cols[i]], flush = True)
            else:
                break
        return data.tolist()
    elif file_form == "tsv":
        df = pd.read_csv(file_name, sep='\t', header = None)
        return list(df[0])
    elif file_form == "gz":
        with gzip.open(file_name, 'rb') as f_in:
            with open(label, 'wb') as f_out:
                shutil.copyfileobj(f_in,f_out)
        read_files(label)
    elif file_form == "csv":
        return pd.read_csv(file_name, sep = ',')
    else:
        print("Sorry the file format is not supported currently, please contact the maintainer if you still want to use it.", file=sys.stderr)

def export_data(gn, Paths):
    os.chdir('./tmp_datasets')
    print(os.getcwd(),flush = True)
    
    print('Identifying file type...',flush = True)
    File_names = []
    for i in Paths:
        File_names.append(i.split('/')[-1])
    print(File_names, flush = True)
    print('Converting to assay file...',flush = True)
    if len(File_names) == 3:
        data = read_files(File_names[0])
        #B = Matrix.todense()
        #df = pd.SparseDataFrame(Matrix)
        #data = df.values.tolist()

        barcodes = read_files(File_names[1])[0:10000]
        genes = read_files(File_names[2])

        exported_assay = {
            "matrix":  data,
            "sampleIds": barcodes,
            "geneIds": genes
        }
    elif len(File_names) == 1:
        pass
    else:
        print("The input files cannot be processed.", file=sys.stderr)
    #print(len(exported_assay["matrix"]), flush=True)
    #print(len(ds.ca["CellID"].tolist()), flush=True)
    gn.export(exported_assay,  "HTAN assay")
    gn.add_result("Successfully exporting HTAN data", data_type='markdown')  
    gn.commit()




def main():
    gn = granatum_sdk.Granatum()

    SIDs = gn.get_arg('SIDs')
    NUM = gn.get_arg("NUM")
    Paths = download_data(SIDs, NUM)
    export_data(gn, Paths)


if __name__ == "__main__":
    main()
