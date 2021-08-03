from deepimpute.multinet import MultiNet
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import sparse
import psutil
import os
import tracemalloc
import csv
from tqdm import tqdm
from time import time
import gc
import json
import gzip
from io import StringIO, BytesIO
import base64


def get_ava_memory():
    usage = psutil.virtual_memory()
    ava_mb = usage.available/(1024*1024)
    print("Ava memory (MB):", ava_mb, flush=True)


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


def run_deep_impute(assay, file):
    # Simulate DeepImpute Gbox
    print("============> Start Simulation  ==============>", flush=True)
    output = {}
    assay = decode_json(assay[file])
    data = np.array(assay.get("matrix")).T
    np.random.seed(12345)
    frameddata = pd.DataFrame(data)
    model = MultiNet()
    model.fit(frameddata, NN_lim="auto", cell_subset=1, minVMR=0.5)
    imputed = model.predict(frameddata, imputed_only=False, policy="restore")
    LABELS_PARAMS = {"fontsize": 14, "fontweight": "bold", "fontname": "serif"}
    data = data + 1
    imputed = imputed + 1
    vmax = np.percentile(np.log10(data.flatten()), 99)
    fig, ax = plt.subplots(1, 2)
    ax[0].imshow(np.log10(data), aspect="auto", vmax=vmax)
    ax[1].imshow(np.log10(imputed), aspect="auto", vmax=vmax)
    ax[0].set_xlabel("Genes", **LABELS_PARAMS)
    ax[1].set_xlabel("Genes", **LABELS_PARAMS)
    ax[0].set_ylabel("Cells", **LABELS_PARAMS)
    ax[0].set_title("raw (log)", **LABELS_PARAMS)
    ax[1].set_title("imputed (log)", **LABELS_PARAMS)
    nb_genes = len(set(model.targets.flatten()))
    def calc_dropout(matrix):
        return np.sum((np.array(matrix) == 0)) * 1. / data.size
    r, p = model.score(frameddata)
    rows, cols = frameddata.shape
    data = data - 1
    imputed = imputed - 1
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
    del data
    gc.collect()
    assay["matrix"] = imputed.T.to_numpy().tolist()
    output["chunk1"] = compress_assay(assay)
    del imputed, assay
    print("============> End Simulation <===========", flush=True)

def main():

    # read data from ./datasets and perform deepimpute
    # save the result to rows
    with open("performance1.csv", "w", buffering=1) as f:
        with open("time_performance1.csv", "w", buffering=1) as f2:
            f.write("Gene Size,Cell Size,Percent,Peak Memory Usage\n")
            f2.write("Gene Size,Cell Size,Percent,Time\n")
            count = 0
            for file in tqdm(sorted(os.listdir("./datasets"))):
                if 0 <= count <50:
                    print(file,flush=True)
                    get_ava_memory()
                    file_path = os.path.join("./datasets", file)
                    
                    assay = json.load(open(file_path, "r"))
                    begin = time()
                    tracemalloc.start()
                    
                    run_deep_impute(assay, file)

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
                count += 1
    f.close()
    f2.close()


if __name__ == "__main__":
    main()
