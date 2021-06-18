import granatum_sdk
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from deepimpute.deepImpute import deepImpute 
import requests
import os
from os.path import basename
import loompy
from tqdm import tqdm
import shutil
import copy
from collections import defaultdict

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



def download_data(ProjectID, gn):
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
                    download_file(url, dest_path)
                    file_urls.add(url)
        print("Finished downloading!", flush = True)
    except requests.exceptions.HTTPError:
        print("Invalid ID entered. Please try again.")
        gn.commit()



def main():
    gn = granatum_sdk.Granatum()

    download_data('4a95101c-9ffc-4f30-a809-f04518a23803', gn)

    os.chdir('./tmp_datasets')
    #print(len(assay['matrix']),flush = True)
    #print(len(assay['sampleIds']), flush = True)

    filename = 't-cell-activation-human-blood-10XV2.loom'
    ds = loompy.connect(filename)

    print('loading assay...', flush = True)

    assay = gn.get_import("assay")
    assay['geneIds'] = ds.ra["Gene"].tolist()
    assay['sampleIds'] = ds.ca["CellID"].tolist()
    data = ds[:,:].T

    seed = gn.get_arg("seed")
    checkbox = gn.get_arg("use_auto_limit")
    cell_subset = gn.get_arg("cell_subset")

    NN_lim = {False: gn.get_arg("NN_lim"), True: "auto"}.get(checkbox, True)

    #model = MultiNet(n_cores="all", seed=seed)
    #model.fit(data, NN_lim=NN_lim, cell_subset=cell_subset)

    frameddata = pd.DataFrame(data)
    imputed, model = deepImpute(frameddata, NN_lim=NN_lim, cell_subset=cell_subset)

    LABELS_PARAMS = {"fontsize": 14, "fontweight": "bold", "fontname": "serif"}

    vmax = np.percentile(np.log10(1 + data.flatten()), 99)

    print("Generating Heatmap")
    fig, ax = plt.subplots(1, 2)
    ax[0].imshow(np.log10(1 + data), aspect="auto", vmax=vmax)
    ax[1].imshow(np.log10(1 + imputed), aspect="auto", vmax=vmax)
    ax[0].set_xlabel("Genes", **LABELS_PARAMS)
    ax[1].set_xlabel("Genes", **LABELS_PARAMS)
    ax[0].set_ylabel("Cells", **LABELS_PARAMS)
    ax[0].set_title("raw (log)", **LABELS_PARAMS)
    ax[1].set_title("imputed (log)", **LABELS_PARAMS)

    gn.add_current_figure_to_results("Heatmaps")
    nb_genes = len(set(model.targets.flatten()))
    #nb_genes = np.sum([len(net.targetGenes) for net in model.targets])

    def calc_dropout(matrix):
        return np.sum((np.array(matrix) == 0)) * 1. / data.size

    r, p = model.score(frameddata)
    rows, cols = frameddata.shape

    message = "\n".join(
        [
            "  - Data frame number of rows: **{0}**",
            "  - Data frame number of columns: **{1}**",
            "  - Number of imputed genes: **{2}**",
            "  - Percentage of dropout entries *before* imputation: **{3:.2f}%**",
            "  - Percentage of dropout entries *after* imputation: **{4:.2f}%**",
            "  - Accuracy (correlation) on masked data: **{5:.2f}**"
        ]
    ).format(
        rows,
        cols,
        nb_genes,
        calc_dropout(data) * 100,
        calc_dropout(imputed.to_numpy()) * 100,
	r
    )

    gn.add_result(message, data_type="markdown")

    print("storing assay...", flush = True)

    assay["matrix"] = imputed.T.to_numpy().tolist()
    gn.export_statically(assay, "Imputed assay")

    gn.commit()


if __name__ == "__main__":
    main()
