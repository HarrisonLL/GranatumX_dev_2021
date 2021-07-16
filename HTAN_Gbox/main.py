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
import gc
import json
import zipfile

# Function to compress chunks
def zipDir(dirpath, outFullName):
    zip = zipfile.ZipFile(outFullName, "w", zipfile.ZIP_DEFLATED)
    for path, dirnames, filenames in os.walk(dirpath):
        fpath = path.replace(dirpath, '')

        for filename in filenames:
            zip.write(os.path.join(path, filename), os.path.join(fpath, filename))
    zip.close

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
    # Store file format and name
    file_form = file_name.split(".")[-1]
    label = '.'.join(file_name.split(".")[:-1])
    
    if file_form == "mtx":
        Matrix = (scipy.io.mmread(file_name))
        gene_nums = Matrix.shape[0]
        cell_nums = Matrix.shape[1]
        cols = Matrix.col.tolist()
        rows = Matrix.row.tolist()
        values = Matrix.data.tolist()
        # return coo sparse matrix info
        return gene_nums, cell_nums, values, rows, cols
    elif file_form == "tsv":
        df = pd.read_csv(file_name, sep='\t', header = None)
        # return the first list of the tsv file
        return list(df[0])
    elif file_form == "gz":
        with gzip.open(file_name, 'rb') as f_in:
            with open(label, 'wb') as f_out:
                shutil.copyfileobj(f_in,f_out)
        return read_files(label)
    elif file_form == "csv":
        return pd.read_csv(file_name, sep = ',')
    else:
        print("Sorry the file format is not supported currently, please contact the maintainer if you still want to use it.", file=sys.stderr)

def export_data(gn, Paths):
    os.chdir('./tmp_datasets')
    
    print('Identifying file type...',flush = True)
    File_names = []
    for i in Paths:
        File_names.append(i.split('/')[-1])
    
    print('Converting to assay file...',flush = True)
    if len(File_names) == 3:
        gene_nums, cell_nums, values, rows, cols = read_files(File_names[0])
        #B = Matrix.todense()
        #df = pd.SparseDataFrame(Matrix)
        #data = df.values.tolist()
        start = 0
        step = 5000
        count = 1
        #chunks = []
        chunk_dir = os.path.join(gn.exports_dir,"chunks")
        #print(chunk_dir, flush = True)
        os.mkdir(chunk_dir)
        print(gn.exports_dir, flush = True)
        while start < cell_nums - 1:
            print("Start to export "+str(count)+" chunk", flush = True)
            end = start + step - 1
            if end < cell_nums:
                data = np.zeros((gene_nums,step))
            else:
                end = cell_nums - 1
                data = np.zeros((gene_nums, end - start + 1))
            
            for i in range(len(cols)):
                if cols[i] >= start and cols[i] <= end:
                    data[rows[i]][cols[i]] = values[i]
                else:
                    break

            barcodes = read_files(File_names[1])[start:end + 1]
            genes = read_files(File_names[2])

            exported_assay = {
                "matrix":  data.tolist(),
                "sampleIds": barcodes,
                "geneIds": genes
            }
            print(count, flush = True)
            assay_name = 'HTAN assay' + str(count)
            #print(os.path.join(chunk_dir, assay_name), flush = True)
            with open(os.path.join(chunk_dir, assay_name), "w") as f:
                json.dump(exported_assay, f)
            #chunks.append(json.dumps(exported_assay))
            #gn.export(exported_assay, assay_name,"assay")
            files = os.listdir(chunk_dir)
            print(files, flush = True)

            count += 1
            del data, exported_assay
            gc.collect()
            start = end
            print("Finished!", flush = True)
            #if count == 3:
            #    break
        #gn.export(chunks, "HTAN chunk", "assay")
        output_path = os.path.join(gn.exports_dir, "chunks.zip")
        zipDir(chunk_dir, output_path)
        files = os.listdir(gn.exports_dir)
        print(files, flush = True)
        gn.dynamic_exports.append({"extractFrom": "chunks.zip", "kind": "assay", "meta": None})

    elif len(File_names) == 1:
        pass
    else:
        print("The input files cannot be processed.", file=sys.stderr)
    #print(len(exported_assay["matrix"]), flush=True)
    #print(len(ds.ca["CellID"].tolist()), flush=True)
    #gn.export_statically(exported_assay,  "HTAN assay")
    gn.add_result("Successfully exporting HTAN data", data_type='markdown')
    #print(gn.exports_anno_file, flush = True)
    #print(gn.dynamic_exports, flush = True)
    #gn.dynamic_exports = []
    #gn.dynamic_exports.append({"extractFrom": "chunks", "kind": "assay", "meta": None}) 
    gn.commit()




def main():
    gn = granatum_sdk.Granatum()

    SIDs = gn.get_arg('SIDs')
    NUM = gn.get_arg("NUM")
    Paths = download_data(SIDs, NUM)
    export_data(gn, Paths)


if __name__ == "__main__":
    main()
