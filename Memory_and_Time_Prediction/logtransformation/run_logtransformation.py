#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time

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

def get_ava_memory():
  usage = psutil.virtual_memory()
  ava_mb = usage.available/(1024*1024)
  print("Ava memory (MB):", ava_mb, flush=True)

def run_logtrans(assay):

  # Simulate Normalization Gbox
  print("============> Start Simulation  ==============>", flush=True)
  output = {}
  log_base = 1
  pseudo_counts = 2
  assay = decode_json(assay['chunk1'])
  print(assay.keys(), flush=True)
  matrix = np.array(assay.get("matrix"))
  #print(matrix[0], flush=True)
  transformed_matrix = np.log(matrix + pseudo_counts) / np.log(log_base)
  non_zero_values_before = matrix.flatten()
  non_zero_values_before = non_zero_values_before[(non_zero_values_before > np.percentile(non_zero_values_before, 5))]
  non_zero_values_after = transformed_matrix.flatten()
  non_zero_values_after = non_zero_values_after[(non_zero_values_after > np.percentile(non_zero_values_after, 5))]

  # plt.figure()
  # plt.subplot(2, 1, 1)
  # plt.title('Before log transformation')
  # plt.hist(non_zero_values_before, bins=100)
  # plt.ylabel('Frequency')
  # plt.xlabel('Expression level')

  # plt.subplot(2, 1, 2)
  # plt.title('After log transformation')
  # plt.hist(non_zero_values_after, bins=100)
  # plt.ylabel('Frequency')
  # plt.xlabel('Expression level')

  # plt.tight_layout()
  # print("Finished plotting", flush = True)
  # caption = (
  #     'The distribution of expression level before and after log transformation. Only the values greater '
  #     'than the 5 percentile (usually zero in single-cell data) and lower than 95 percentile are considered.'
  # )

  del matrix, non_zero_values_before, non_zero_values_after
  gc.collect()
  assay["matrix"] = transformed_matrix.tolist()
  output["chunk1"] = compress_assay(assay)
  del transformed_matrix,assay
  gc.collect()
  print("============> End Simulation <===========", flush=True)

def main():

    # read data from ./datasets and perform deepimpute
    # save the result to rows
    count = 0
    with open("mem_performance_logtrans.csv", "a", buffering=1) as f:
      with open("time_performance2_logtrans.csv", "a", buffering=1) as f2:

        #f.write("Gene Size,Cell Size,Percent,Peak Memory Usage\n")
        #f2.write("Gene Size,Cell Size,Percent,Time\n")
        for file in tqdm(sorted(os.listdir("/content/drive/MyDrive/datasets"))):
          print(file,flush=True)
          count += 1
          while count > 56:
            get_ava_memory()
            file_path = os.path.join("/content/drive/MyDrive/datasets", file)
            
            assay = json.load(open(file_path, "r"))
            begin = time()
            tracemalloc.start()
            
            run_logtrans(assay)

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