from itertools import combinations
import scipy
import multiprocessing
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
from sklearn.preprocessing import quantile_transform
from scipy.sparse import csc_matrix
import time
from granatum_sdk import Granatum
from granatum_sdk_subclass import granatum_extended
from granatum_sdk_subclass2 import granatum_extended2
# import pandas as pd
# import seaborn as sns


nans = np.array([np.nan, np.nan])
zeros = np.array([0, 0])


def trim_extreme(x, a, b):
    low = np.percentile(x, a)
    high = np.percentile(x, b)
    filtered = x[(x > low) & (x < high)]
    return filtered.copy()


def make_plot(adata, log_trans=False):
    violin_data = []
    for cell in adata.X:
        filtered = cell.toarray().flatten()
        #filtered = trim_extreme(filtered, 5, 95)
        if log_trans:
            #cell = np.log1p(cell)
            filtered = np.log1p(filtered)
        if filtered.shape[0] == 0:
            #cell = zeros
            filtered = zeros

        violin_data.append(filtered)

    plt.figure()
    plt.boxplot(violin_data)
    plt.xlabel('Cells')
    plt.ylabel('Expression lvl (log transformed)')
    plt.tight_layout()

def quantile_normalization(mat):
    # double argsort for getting the corresponding ranks for
    # each element in the vector

    rank_mat = np.argsort(np.argsort(mat, 1), 1)
    medians = np.median(np.sort(mat, 1), 0)
    normalized = np.zeros_like(mat)

    for i in range(rank_mat.shape[0]):
       normalized[i, :] = medians[rank_mat[i, :]]

    # normalized = quantile_transform(mat, copy=False)

    #return normalized.tolist()
    return sc.AnnData(csc_matrix(normalized))


def main():
    bool_sparse = False
    gn = granatum_extended("scanpynormalization")
    chunks = gn.get_import('assay')
    if chunks["current chunk"][-1] == "sparse":
        gn = granatum_extended2("scanpygenefilering")
        chunks = gn.get_import('assay')
        bool_sparse = True
    output = {"origin data size":chunks["origin data size"], "current chunk":["scanpynormalization", "row", "sparse"], "suggested chunk":chunks["suggested chunk"]}
    gn.adjust_transform(chunks)
    switch_time = time.time()
    print("--- %s seconds --- for transforming" % (switch_time - start_time), flush = True)
    chunks = gn.new_file

    for i in range(len(chunks)):
        combined = gn.combine_new_chunk(chunks["chunk"+str(i+1)], "col")
        if bool_sparse:
            matrix = scipy.sparse.csc_matrix((combined.get("data"), combined.get("indices"), combined.get("indptr")), shape=(len(combined["geneIds"]), len(combined["sampleIds"]))).todense()
            chunk1_assay["matrix"] = matrix
            del(combined["data"])
            del(combined["indices"])
            del(combined["indptr"])

    #chunk1 = gn.decompress_chunk(chunks["chunk1"])
    #matrix = chunk1['matrix']
    #matrix = scipy.sparse.csc_matrix((chunk1.get("data"), chunk1.get("indices"), chunk1.get("indptr")), shape=(len(chunk1["geneIds"]), len(chunk1["sampleIds"]))).todense()

        adata = gn.ann_data_from_assay(combined)
        num_cells_to_sample = gn.get_arg('num_cells_to_sample')
        method = gn.get_arg('method')
        log_trans_when_plot = gn.get_arg('log_trans_when_plot')

        if num_cells_to_sample > adata.shape[0]:
            num_cells_to_sample = adata.shape[0]

        sampled_cells_idxs = np.sort(np.random.choice(adata.shape[0], num_cells_to_sample, replace=False))
        if i == len(chunks):
            make_plot(adata[sampled_cells_idxs, :], log_trans=log_trans_when_plot)
            gn.add_current_figure_to_results(
                'Before normalization: Each bar in the box plot represents one cell.',
                height=350,
                dpi=75 * 40 / max(40, num_cells_to_sample)
            )

        if method == 'quantile':
            adata2 = quantile_normalization(adata.X.toarray())
            adata2.var_names = adata.var_names.tolist()
            adata2.obs_names = adata.obs_names.tolist()
            adata = adata2
        elif method == 'scanpy':
            sc.pp.normalize_total(adata)
        else:
            raise ValueError()
        if i == len(chunks):
            make_plot(adata[sampled_cells_idxs, :], log_trans=log_trans_when_plot)
            gn.add_current_figure_to_results(
                'After normalization: Each bar in the box plot represents one cell.',
                height=350,
                dpi=75 * 40 / max(40, num_cells_to_sample)
            )
        print("Compressing", flush = True)
        output["chunk"+str(i+1)] = gn.compress_chunk(gn.assay_from_ann_data(adata))
       
    gn.export(output, 'Normalized assay',dynamic=False)
    gn.commit()
    print("--- %s seconds --- for whole gbox" % (time.time() - start_time), flush = True)
    time.sleep(20)


if __name__ == '__main__':
    main()
