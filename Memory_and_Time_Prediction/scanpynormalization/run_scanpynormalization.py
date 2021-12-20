#!/usr/bin/env python

from itertools import combinations
import multiprocessing
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
from sklearn.preprocessing import quantile_transform
from scipy.sparse import coo_matrix
from scipy.sparse import csc_matrix

from tqdm import tqdm
from time import time
import gc
import json
import gzip
from io import StringIO, BytesIO
import base64
import psutil
import os
import tracemalloc
# import pandas as pd
# import seaborn as sns


nans = np.array([np.nan, np.nan])
zeros = np.array([0, 0])

def decode_json(jsonstr):
    bio = BytesIO()
    stream = BytesIO(base64.b64decode(jsonstr.encode('utf-8')))
    decompressor = gzip.GzipFile(fileobj=stream, mode='r')

    while True:
        chunk = decompressor.read(8192)
        if not chunk:
            decompressor.close()
            bio.seek(0)
            output = bio.read().decode("utf-8")
            break
        bio.write(chunk)
    return json.loads(output)


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
    plt.close()

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

def get_ava_memory():
  usage = psutil.virtual_memory()
  ava_mb = usage.available/(1024*1024)
  print("Ava memory (MB):", ava_mb, flush=True)

def run_scnorm(assay, file):

  # Simulate Normalization Gbox
  print("============> Start Simulation  ==============>", flush=True)
  output = {}
  assay = decode_json(assay['chunk1'])
  sparse_matrix = coo_matrix(assay.get("matrix")).tocsc()
  adata = sc.AnnData(sparse_matrix.transpose())

  adata.var_names = assay.get("geneIds")
  adata.obs_names = assay.get("sampleIds")

  num_cells_to_sample = 40
  method = 'quantile'
  log_trans_when_plot = 0

  if num_cells_to_sample > adata.shape[0]:
    num_cells_to_sample = adata.shape[0]

  sampled_cells_idxs = np.sort(np.random.choice(adata.shape[0], num_cells_to_sample, replace=False))
  make_plot(adata[sampled_cells_idxs, :], log_trans=log_trans_when_plot)

  if method == 'quantile':
    adata2 = quantile_normalization(adata.X.toarray())
    adata2.var_names = adata.var_names.tolist()
    adata2.obs_names = adata.obs_names.tolist()
    adata = adata2
  elif method == 'scanpy':
    sc.pp.normalize_total(adata)
  else:
    raise ValueError()

  make_plot(adata[sampled_cells_idxs, :], log_trans=log_trans_when_plot)

  del adata2
  gc.collect()
  assay["matrix"] = adata.X.T.toarray().tolist()
  output["chunk1"] = compress_assay(assay)
  del adata, assay
  gc.collect()
  print("============> End Simulation <===========", flush=True)

def main():

    # read data from ./datasets and perform deepimpute
    # save the result to rows
    with open("performance_scnorm.csv", "w", buffering=1) as f:
      with open("time_performance_scnorm.csv", "w", buffering=1) as f2:
        f.write("Gene Size,Cell Size,Percent,Peak Memory Usage\n")
        f2.write("Gene Size,Cell Size,Percent,Time\n")
        for file in tqdm(sorted(os.listdir("../../../GranatumX_dense_datasets"))):
          print(file,flush=True)
          get_ava_memory()
          file_path = os.path.join("../../../GranatumX_dense_datasets", file)
          
          assay = json.load(open(file_path, "r"))
          begin = time()
          tracemalloc.start()
          
          run_scnorm(assay, file)

          _,PEAK = tracemalloc.get_traced_memory()
          tracemalloc.stop()
          end = time()
          elapse = end - begin


          print("Peak Usage of the function:" + str(PEAK/1024/1024),flush=True)
          tmp = file.replace(".gz", "").split("_")
          gene_size = tmp[0]
          cell_size = tmp[1]
          percent = tmp[2]
          print("Finished simualation for sample:", flush=True)
          print(gene_size, cell_size, percent, flush=True)
          f.write(",".join([gene_size, cell_size, percent, str(PEAK/1024/1024)])+ "\n")
          f2.write(",".join([gene_size, cell_size, percent, str(elapse)])+ "\n")
          del(PEAK)
          del(tmp)
          del(gene_size)
          del(cell_size)
          del(percent)
          del(begin)
          del(end)
          del(elapse)
          gc.collect()
    f.close()
    f2.close()


if __name__ == "__main__":
    main()
