import granatum_sdk
import os
from os.path import basename
import loompy
from tqdm import tqdm
import shutil
import copy
import synapseclient
import pandas as pd
import numpy as np


def download_file(url, output_path):
    url = url.replace('/fetch', '')  # Work around https://github.com/DataBiosphere/azul/issues/2908
    
    response = requests.get(url, stream=True)
    response.raise_for_status()
    
    total = int(response.headers.get('content-length', 0))
    print(f'Downloading to: {output_path}', flush=True)
    
    with open(output_path, 'wb') as f:
        with tqdm(total=total, unit='B', unit_scale=True, unit_divisor=1024) as bar:
            for chunk in response.iter_content(chunk_size=1024):
                size = f.write(chunk)
                bar.update(size)


def export_data(gn):
    os.chdir('./tmp_datasets')
    #print(len(assay['matrix']),flush = True)
    #print(len(assay['sampleIds']), flush = True)
    ds = None
    count = 0
    for file in os.listdir("./"):
        if file.endswith(".loom"):
            ds = loompy.connect(file)
            count += 1
        if count == 2:
            print(file)
            break
            
    print('Converting to assay file...',flush = True)
    length = len(ds[0, :])
    print('Choose 1000 cells out of %i'%length, flush=True)
    indices = sorted(np.random.choice(length,1000, replace=False))
    exported_assay = {
        "matrix":  (ds[:,indices].tolist()),
        "sampleIds": ((ds.ca["CellID"])[indices].tolist()),
        "geneIds": ds.ra["Gene"].tolist(),
    }
    #print(len(exported_assay["matrix"]), flush=True)
    #print(len(ds.ca["CellID"].tolist()), flush=True)
    gn.export(exported_assay,  "HCA assay")
    gn.add_result("Successfully exporting HCA data", data_type='markdown')  
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
    # TO-DO:
    # modify the function to make export correctly
    # export_data(gn)
    gn.commit()


def main():
    gn = granatum_sdk.Granatum()

    SID = gn.get_arg('SID')
    download_data(SID, gn)


if __name__ == "__main__":
    main()
