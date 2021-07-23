import granatum_sdk
import requests
import os
from os.path import basename
import loompy
from tqdm import tqdm
import shutil
import copy
import pandas as pd
import numpy as np
import scipy.io
import scipy.sparse
from sympy import symbols, Eq, solve
from time import time
import psutil
import zipfile
import json
import gc
from io import StringIO, BytesIO
import gzip
import base64


def pred_cell_size(coef, genesize, percent, ava_mem):
    x = symbols('x')
    eq = Eq(coef[0] + coef[1]*genesize + coef[2]*x + coef[3]*percent + coef[4]*(genesize**2) + coef[5]*(genesize*x) + coef[6]*(genesize*percent)+ coef[7]*(x**2) + coef[8]*(x*percent)+coef[9]*(percent**2)-ava_mem)
    print(solve(eq), flush=True)
    return [int(sol) for sol in solve(eq) if sol > 0]


# Function to compress chunks
def compress_assay(exported_assay):
    tmp = json.dumps(exported_assay)
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


def iterate_matrices_tree(tree, keys=()):
    if isinstance(tree, dict):
        for k, v in tree.items():
            yield from iterate_matrices_tree(v, keys=(*keys, k))
    elif isinstance(tree, list):
        for file in tree:
            yield keys, file
    else:
        assert False


def export_data(gn,coeffs,ava_mem):
    os.chdir('./tmp_datasets')
    ds = loompy.connect(os.listdir("./")[0])
    cell_size = ds.shape[1]
    gene_size = ds.shape[0]
    sparse = ds.sparse()
    percent = np.sum(sparse.data < 101) * 100 / (ds.shape[0] * ds.shape[1])
    print(ava_mem, flush=True)
    print(cell_size, flush=True)
    print(gene_size, flush=True)
    chunk_size = 1000
    print(chunk_size,flush=True)

    if chunk_size >= cell_size:
        print("Exporting all data..", flush=True)
        # TO-DO
        exported_assay = {
                "matrix":  (ds[:,:].tolist()),
                "sampleIds": ((ds.ca["CellID"])[indices].tolist()),
                "geneIds": ds.ra["Gene"].tolist(),
        }
        gn.export(exported_assay,  "HCA assay")
    else:
        print("Chunking the data..", flush=True)
        output = {}
        for i in range(0, cell_size, chunk_size):
            start = time()
            data = sparse.tocsr()[:,i:i+chunk_size].todense().tolist()
            end = time()
            print(end-start, flush=True)
            exported_assay = {
                "matrix":  data,
                "sampleIds": ((ds.ca["CellID"])[i:i+chunk_size].tolist()),
                "geneIds": ds.ra["Gene"].tolist(),
            }
            chunk_dir = os.path.join(gn.exports_dir,"chunks")
            assay_name = "HCA assay" + str(i//10) + '.gz'
            output[assay_name] = compress_assay(exported_assay)

            del(data)
            del(exported_assay)
            gc.collect()
        
        with gzip.open(os.path.join(gn.exports_dir,"chunks.gz"),"wt") as f:
            json.dump(output, f)
        gn.dynamic_exports.append({"extractFrom": "chunks.gz", "kind": "assay", "meta": None})

   #print('Converting to assay file...',flush = True)
   # length = len(ds[0, :])
   # print('Choose 1000 cells out of %i'%length, flush=True)
   # indices = sorted(np.random.choice(length,1000, replace=False))
   # exported_assay = {
   #     "matrix":  (ds[:,indices].tolist()),
   #     "sampleIds": ((ds.ca["CellID"])[indices].tolist()),
   #     "geneIds": ds.ra["Gene"].tolist(),
   # }
   # gn.export(exported_assay,  "HCA assay")
    gn.add_result("Successfully exporting HCA data", data_type='markdown')  
    gn.commit()


def download_data(ProjectID, Species, Organ, gn,coeffs,ava_mem):
    # destroy if exits, and then create ./tmp_datasets directory
    dirpath = './tmp_datasets'
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
    os.mkdir(dirpath)

    project_uuid = ProjectID
    catalog = 'dcp6'
    endpoint_url = f'https://service.azul.data.humancellatlas.org/index/projects/{project_uuid}'
    save_location = dirpath


    try:
        response = requests.get(endpoint_url, params={'catalog': catalog})
        response.raise_for_status()

        response_json = response.json()
        project = response_json['projects'][0]

        file_urls = set()
        for key in ('matrices', 'contributorMatrices'):
            tree = project[key]
            for path, file_info in iterate_matrices_tree(tree):
                url = file_info['url']
                if url not in file_urls:
                    dest_path = os.path.join(save_location, file_info['name'])
                    # only download DCP2-processed files with selected organ type and species
                    name = file_info['name']
                    if  (name.endswith('.loom')) and Organ in name.lower() and Species in name.lower():
                        download_file(url, dest_path)
                    file_urls.add(url)
        
        directory= os.listdir('./tmp_datasets')
        if len(directory) == 0:
            print("Download unsucessful. The file of interest is not found in the project folder.", flush=True)
            print("Please check if any mistyping occurs.", flush=True)
            gn.commit()
            return
        
        print("Finished downloading!", flush = True)
        export_data(gn,coeffs,ava_mem)
    
    except requests.exceptions.HTTPError:
        print("Invalid ID entered. Please try again.")
        gn.commit()


def main():
    df = pd.read_csv("coeffs.csv")
    coeffs = df.iloc[:1, 1:].values.squeeze().tolist()
    print(coeffs, flush=True)
    ava_mem = (psutil.virtual_memory()).available/(1024*1024)

    gn = granatum_sdk.Granatum()
    ProjectID = gn.get_arg('PID').lower()
    Species = gn.get_arg('SPE').lower()
    Organ = gn.get_arg('ORG').lower()
    download_data(ProjectID, Species, Organ, gn,coeffs,ava_mem)


if __name__ == "__main__":
    main()
