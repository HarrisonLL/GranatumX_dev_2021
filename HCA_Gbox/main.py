import granatum_sdk
import requests
import os
from os.path import basename
import loompy
from tqdm import tqdm
import shutil
import copy
from collections import defaultdict
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


def iterate_matrices_tree(tree, keys=()):
    if isinstance(tree, dict):
        for k, v in tree.items():
            yield from iterate_matrices_tree(v, keys=(*keys, k))
    elif isinstance(tree, list):
        for file in tree:
            yield keys, file
    else:
        assert False


def export_data(gn):
    os.chdir('./tmp_datasets')
    ds = loompy.connect(os.listdir("./")[0])
    print('Converting to assay file...',flush = True)
    length = len(ds[0, :])
    print('Choose 1000 cells out of %i'%length, flush=True)
    indices = sorted(np.random.choice(length,1000, replace=False))
    exported_assay = {
        "matrix":  (ds[:,indices].tolist()),
        "sampleIds": ((ds.ca["CellID"])[indices].tolist()),
        "geneIds": ds.ra["Gene"].tolist(),
    }
    gn.export(exported_assay,  "HCA assay")
    gn.add_result("Successfully exporting HCA data", data_type='markdown')  
    gn.commit()


def download_data(ProjectID, Species, Organ, gn):
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
        export_data(gn)
    
    except requests.exceptions.HTTPError:
        print("Invalid ID entered. Please try again.")
        gn.commit()


def main():
    gn = granatum_sdk.Granatum()
    ProjectID = gn.get_arg('PID').lower()
    Species = gn.get_arg('SPE').lower()
    Organ = gn.get_arg('ORG').lower()
    download_data(ProjectID, Species, Organ, gn)


if __name__ == "__main__":
    main()
