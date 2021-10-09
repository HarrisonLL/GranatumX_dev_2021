#!/usr/bin/env python

from itertools import combinations

import multiprocessing

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

import scipy
import gc
import time
from granatum_sdk_subclass import granatum_extended
from granatum_sdk_subclass2 import granatum_extended2

# import pandas as pd
# import seaborn as sns


def to_percentage(x):
    return "%{:3.2f}".format(x * 100)

def main():
    
    gn = granatum_extended("scanpypca")
    chunks = gn.get_import('assay')
    org_chunk_kind = chunks["current chunk"][-2]
    threshold = int(chunks["suggested chunk"]["scanpypca"][-1])
    # for testing use
    threshold = 10000

    bool_sparse = False
    if chunks["current chunk"][-1] == "sparse":
        gn = granatum_extended2("scanpypca")
        chunks = gn.get_import('assay')
        bool_sparse = True
    
    # convert back to full sparse dataset
    combined = gn.combine_new_chunk([chunks["chunk"+str(i+1)] for i in range(len(chunks)-3)],org_chunk_kind)
    if bool_sparse:
        if org_chunk_kind == "col":
            sparse_matrix  = scipy.sparse.csc_matrix((combined.get("data"), combined.get("indices"), combined.get("indptr")), shape=(len(combined["geneIds"]), len(combined["sampleIds"])))
        else:
            sparse_matrix  = scipy.sparse.csr_matrix((combined.get("data"), combined.get("indices"), combined.get("indptr")), shape=(len(combined["geneIds"]), len(combined["sampleIds"]))).tocsc()
        adata = sc.AnnData(sparse_matrix.transpose())
        adata.var_names = combined.get("geneIds")
        adata.obs_names = combined.get("sampleIds")

    else:
        adata = gn.ann_data_from_assay(combined)
    del combined
    gc.collect()
    

    print(adata.shape, flush=True)
    num_top_comps = gn.get_arg("num_top_comps")

    if adata.shape[0] < threshold:  
        start1 = time.time()
        print("starting to do PCA", flush=True)
        sc.pp.pca(adata, 20))
    else:
        start1 = time.time()
        print("starting to do IPCA", flush=True)
        sc.pp.pca(adata, 20, chunked=True, chunk_size=1000) # NEED TO TEST CHUNK SIZE
    
    variance_ratios = adata.uns["pca"]["variance_ratio"]
    pc_labels = ["PC{}".format(x + 1) for x in range(len(variance_ratios))]

    plt.figure()
    plt.bar(pc_labels, variance_ratios)
    plt.tight_layout()
    gn.add_current_figure_to_results("Explained variance (ratio) by each Principal Component (PC)", height=350, dpi=75)

    X_pca = adata.obsm["X_pca"]

    for i, j in combinations(range(num_top_comps), 2):
        xlabel = "PC{}".format(i + 1)
        ylabel = "PC{}".format(j + 1)

        plt.figure()
        plt.scatter(X_pca[:, i], X_pca[:, j], s=5000 / adata.shape[0])
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.tight_layout()
        gn.add_current_figure_to_results("PC{} vs. PC{}".format(i + 1, j + 1), dpi=75)

        pca_export = {
            "dimNames": [xlabel, ylabel],
            "coords": {sample_id: X_pca[k, [i, j]].tolist() for k, sample_id in enumerate(adata.obs_names)},
        }
        gn.export(pca_export, "PC{} vs. PC{}".format(i + 1, j + 1), kind="sampleCoords", meta={})
    print(time.time()-start1, flush=True)
    gn.commit()


if __name__ == "__main__":
    main()
