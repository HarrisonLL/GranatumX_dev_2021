from granatum_sdk_subclass import granatum_extended
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from deepimpute.deepImpute import deepImpute 
from deepimpute.multinet import MultiNet
import time
import json
import gzip
from io import StringIO, BytesIO
import base64
import gc
import os

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


def main():
    start_time = time.time()
    gn = granatum_extended("deepimpute")
    #output = {}
    
    #firstly save to chunks
    chunks = gn.get_import("assay")
    print(chunks.keys(), flush = True)
    begin = time.time()
    gn.adjust_transform(chunks)
    end = time.time()
    print("Adjusting from 1000 to 3456 took %.3f" %(end-begin), flush=True)
    chunks = gn.new_file
    
    print(len(chunks.keys()), flush = True)
    print(list(chunks.keys()), flush=True)
    
    seed = gn.get_arg("seed")
    checkbox = gn.get_arg("use_auto_limit")
    cell_subset = gn.get_arg("cell_subset")
    NN_lim = {False: gn.get_arg("NN_lim"), True: "auto"}.get(checkbox, True)
    model = MultiNet(seed=seed)
    output = dict()

    for i in range(len(chunks)):
        combined = gn.combine_new_chunk(chunks["chunk"+str(i+1)])
        data = np.array(combined.get("matrix")).T
        frameddata = pd.DataFrame(data)
        model.fit(frameddata, NN_lim=NN_lim, cell_subset=cell_subset)
        del combined, data, frameddata
        gc.collect()
        if i == 1:
            break
    


    def calc_dropout(matrix):
        return np.sum((np.array(matrix) == 0)) * 1. / data.size
    
    
    nb_genes = 0
    sum_r = 0
    for i in range(len(chunks)):
        combined = gn.combine_new_chunk(chunks["chunk"+str(i+1)])
        data = np.array(combined.get("matrix")).T
        frameddata = pd.DataFrame(data)
        imputed = model.predict(frameddata, imputed_only=False, policy="restore")
    
            
        combined["matrix"] = imputed.T.to_numpy().tolist()
        assay = compress_assay(combined)
        output["chunk" + str(i)] = assay
        nb_genes += len(set(model.targets.flatten()))
        r, _ = model.score(frameddata)
        sum_r += r
        rows, cols = frameddata.shape

        if i == 1:
            cols *= 2
            sum_r /= 2
            fig, ax = plt.subplots(1, 2)
            LABELS_PARAMS = {"fontsize": 14, "fontweight": "bold", "fontname": "serif"}

            message = "\n".join(
           [
            "  - Data frame number of rows: **{0}**",
            "  - Data frame number of columns: **{1}**",
            "  - Number of imputed genes: **{2}**",
            "  - Percentage of dropout entries *before* imputation: **{3:.2f}%**",
            "  - Percentage of dropout entries *after* imputation: **{4:.2f}%**",
            "  - Averaged accuracy (correlation) on masked data: **{5:.2f}**"
           ]
            ).format(
                sum_r,
                cols,
                nb_genes,
                calc_dropout(data) * 100,
                calc_dropout(imputed.to_numpy()) * 100,
                r
            )
            gn.add_result(message, data_type="markdown")
            del combined,data,frameddata,imputed
            gc.collect()
            break
        
        del combined,data,frameddata,imputed
        gc.collect()

    gn.export_statically(output, "Imputed assay")
    gn.commit()

    #print("--- %s seconds --- for compress" % (time.time() - deepimpute_time), flush = True)
    #print("--- %s seconds --- for whole gbox" % (time.time() - start_time), flush = True)
    #time.sleep(10)

if __name__ == "__main__":
    main()
