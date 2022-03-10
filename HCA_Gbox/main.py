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
import time

"""
This function returns a predicted cell size
Predicted chunk size is used for DeepImpute module

Notice current coeffcients will not output valid cell size, IF ava_mem < 50000MB
In the above case, we use default size of 1000
"""
def pred_cell_size_degree2(coef,genesize, percent,ava_mem):
    x = symbols('x',real=True)
    eq1 = Eq(coef[0] + coef[1]*genesize + coef[2]*x + coef[3]*percent + coef[4]*(genesize**2) + coef[5]*(genesize*x) + coef[6]*(genesize*percent)+ coef[7]*(x**2) + coef[8]*(x*percent)+coef[9]*(percent**2)-ava_mem)
    #print(solve(eq1), flush=True)
    if len(solve(eq1)) == 0:
        return 1000 # DEFAULT RETURN, CHANGE IT IF NEEDED
    ret = max(solve(eq1))
    return int(ret) if ret > 0 else 1000 # DEFAULT RETURN, CHANGE IT IF NEEDED


"""
This function returns a predicted gene size
Predicted chunk size is used for Genefilter module
"""
def pred_gene_size_degree2(coef,cellsize, percent,ava_mem):
    x = symbols('x',real=True)
    eq1 = Eq(coef[0] + coef[1]*x + coef[2]*cellsize + coef[3]*percent + coef[4]*(x**2) + coef[5]*(x*cellsize) + coef[6]*(x*percent)+ coef[7]*(cellsize**2) + coef[8]*(cellsize*percent)+coef[9]*(percent**2)-ava_mem)
    if len(solve(eq1)) == 0:
        return 1000 # DEFAULT RETURN, CHANGE IT IF NEEDED
    ret = min(solve(eq1)) if min(solve(eq1)) > 0 else max(solve(eq1))
    return int(ret) if ret > 0 else 1000 # DEFAULT RETURN, CHANGE IT IF NEEDED


"""
This function compresses chunks
"""
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

"""
Dataset chunking, chunk size prediction for next modules
"""
def export_data(gn,coeffs,coeffs2,ava_mem):
    os.chdir('./tmp_datasets')
    ds = loompy.connect(os.listdir("./")[0])
    cell_size = ds.shape[1]
    gene_size = ds.shape[0]
    sparse = ds.sparse()
    percent = np.sum(sparse.data < 101) * 100 / (ds.shape[0] * ds.shape[1])
    print("available memory: {}".format(ava_mem), flush=True)
    print("cell_size: {}".format(cell_size), flush=True)
    print("gene_size: {}".format(gene_size), flush=True)
    chunk_size = 1000
    print("default chunk size: {}".format(chunk_size), flush=True)
    deepimpute_size = pred_cell_size_degree2(coeffs, gene_size, percent, ava_mem)
    print("predicted cell size for DeepImpute: {}".format(deepimpute_size),flush=True)
    genefilter_size = pred_gene_size_degree2(coeffs2, cell_size, percent, ava_mem)
    print("predicted gene size for GeneFilter: {}".format(genefilter_size),flush=True)

    print("Chunking the data..", flush=True)
    
    output = {"origin data size":[gene_size, cell_size],
                  "current chunk":["test", "col", "sparse"],
                  "suggested chunk": {
                                    "test":["col",1000],
                                    "deepimpute": ["col", deepimpute_size],
                                    "genefilter":["row", genefilter_size]
    
                                    }
                  }
    count = 1
    for i in range(0, cell_size, chunk_size):
        data = sparse.tocsc()[:,i:i+chunk_size]
        exported_assay = {        
                "sampleIds":((ds.ca["CellID"])[i:i+chunk_size].tolist()),
                "geneIds": ds.ra["Gene"].tolist(),
                "data":data.data.tolist(),
                "indices":data.indices.tolist(),
                "indptr": data.indptr.tolist()
        }
    
        assay_name = "chunk" + str(count)
        output[assay_name] = exported_assay
        count += 1
        del(data)
        del(exported_assay)
        gc.collect()

    with open(os.path.join(gn.exports_dir,"chunks"),"wt") as f:
        json.dump(output, f)
    gn.dynamic_exports.append({"extractFrom": "chunks", "kind": "assay", "meta": None})

    gn.add_result("Successfully exporting HCA data", data_type='markdown')  
    gn.commit()


"""
Download data to server
"""
def download_data(ProjectID, Species, Organ, gn,coeffs,coeffs2,ava_mem):
    # destroy if exits, and then create ./tmp_datasets directory
    dirpath = './tmp_datasets'
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
    os.mkdir(dirpath)

    project_uuid = ProjectID
    endpoint_url = f'https://service.azul.data.humancellatlas.org/index/projects/{project_uuid}'
    save_location = dirpath


    try:
        response = requests.get(endpoint_url)
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
            time.sleep(30)
            gn.commit()
            return
        
        print("Finished downloading!", flush = True)
        export_data(gn,coeffs,coeffs2,ava_mem)
    
    except requests.exceptions.HTTPError:
        print("Invalid ID entered. Please try again.", flush=True)
        time.sleep(30)
        gn.commit()


def main():
    df = pd.read_csv("coeffs_deepimpute.csv")
    df2 = pd.read_csv("coeffs_genefilter.csv")
    coeffs = df.iloc[:1, 1:].values.squeeze().tolist()
    coeffs2 = df2.iloc[:1, 1:].values.squeeze().tolist()
    ava_mem = (psutil.virtual_memory()).available/(1024*1024)

    gn = granatum_sdk.Granatum()
    ProjectID = gn.get_arg('PID').lower()
    Species = gn.get_arg('SPE')
    Organ = gn.get_arg('ORG').lower()
    if Species == "Mus musculus":
        Species = "mouse"
    elif Species == "Homo sapiens":
        Species = "human"
    download_data(ProjectID, Species, Organ, gn,coeffs,coeffs2,ava_mem)


if __name__ == "__main__":
    main()
