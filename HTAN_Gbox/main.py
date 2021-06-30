import granatum_sdk
import os
#from os.path import basename
#import loompy
from tqdm import tqdm
import shutil
import copy
import synapseclient
import pandas as pd
import numpy as np
import scipy.io
import scipy.sparse


def export_data(gn):
    os.chdir('./tmp_datasets')
    print(os.getcwd(),flush = True)
    #print(len(assay['matrix']),flush = True)
    #print(len(assay['sampleIds']), flush = True)
    #ds = None
    #count = 0
    #for file in os.listdir("./"):
    #    if file.endswith(".loom"):
    #        ds = loompy.connect(file)
    #        count += 1
    #    if count == 2:
    #        print(file)
    #        break
    
    print('Converting to assay file...',flush = True)
    #length = len(ds[0, :])
    Matrix = (scipy.io.mmread('HTAPP-272-SMP-4831_none_channel1_matrix.mtx'))
    #B = Matrix.todense()
    df = pd.SparseDataFrame(Matrix)
    data = df.values.tolist()
    #print('Choose 1000 cells out of %i'%length, flush=True)
    #indices = sorted(np.random.choice(length,1000, replace=False))

    barcodes = pd.read_csv("HTAPP-272-SMP-4831_none_channel1_barcodes.tsv", sep='\t', header = None) 
    genes = pd.read_csv("HTAPP-272-SMP-4831_none_channel1_genes.tsv", sep='\t', header = None)

    exported_assay = {
        "matrix":  data,
        "sampleIds": list(barcodes[0]),
        "geneIds": list(genes[0])
    }
    #print(len(exported_assay["matrix"]), flush=True)
    #print(len(ds.ca["CellID"].tolist()), flush=True)
    gn.export(exported_assay,  "HTAN assay")
    gn.add_result("Successfully exporting HTAN data", data_type='markdown')  
    gn.commit()


def download_data(SID, gn):
    # destroy if exits, and then create ./tmp_datasets directory
    dirpath = './tmp_datasets'
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
    os.mkdir(dirpath)

    syn = synapseclient.Synapse()
    syn.login("HarrisL", "Hlld@0217")

    # TO-DO:
    # add try except here for invalid ID
    entity = syn.get("syn24181451", downloadLocation=dirpath)
    entity2 = syn.get("syn24181449", downloadLocation=dirpath)
    entity3 = syn.get("syn24181474", downloadLocation=dirpath)
    print("Download barcode file path is:" + entity.path, flush=True)
    print("Download gene file path is:" + entity2.path, flush=True) 
    print("Download matrix file path is:" + entity3.path, flush=True) 
    # TO-DO:
    # modify the function to make export correctly
    export_data(gn)
    gn.commit()


def main():
    gn = granatum_sdk.Granatum()

    SID = gn.get_arg('SID')
    download_data(SID, gn)


if __name__ == "__main__":
    main()
