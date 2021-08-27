from granatum_sdk_subclass2 import granatum_extended2
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
import scipy

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

    gn = granatum_extended2("deepimpute") 
    #firstly save to chunks
    chunks = gn.get_import("assay")
    print(chunks.keys(), flush = True)
    # create an output
    output = {"origin data size":chunks["origin data size"], "current chunk":["deepimpute", "col"], "suggested chunk":chunks["suggested chunk"]}
    begin = time.time()
    gn.adjust_transform(chunks)
    end = time.time()
    print("Adjusting from 1000 to 500 took %.3f" %(end-begin), flush=True)
    chunks = gn.new_file
    
    print(len(chunks.keys()), flush = True)
    print(list(chunks.keys()), flush=True)
    
    seed = gn.get_arg("seed")
    checkbox = gn.get_arg("use_auto_limit")
    cell_subset = gn.get_arg("cell_subset")
    NN_lim = {False: gn.get_arg("NN_lim"), True: "auto"}.get(checkbox, True)
    model = MultiNet(seed=seed)



    nb_genes = 0
    sum_r = 0
    data_dropout = 0
    impu_dropout = 0
    def calc_dropout(matrix):
        return np.sum((np.array(matrix) == 0)) * 1. / data.size
    
    for i in range(len(chunks)):
        combined = gn.combine_new_chunk(chunks["chunk"+str(i+1)], "col")
        matrix  = scipy.sparse.csc_matrix((combined.get("data"), combined.get("indices"), combined.get("indptr")), shape=(len(combined["geneIds"]), len(combined["sampleIds"]))).todense()
        data = matrix.T
        print(data.shape,flush=True)
        del matrix
        gc.collect()
        
        frameddata = pd.DataFrame(data)
        model.fit(frameddata, NN_lim=NN_lim, cell_subset=cell_subset)
        imputed = model.predict(frameddata, imputed_only=False, policy="restore")
    
        #if i == 0:
        #    tmp = scipy.sparse.csc_matrix(imputed.T.to_numpy())
        #    assay = {
        #        "data":tmp.data,
        #        "indices":tmp.indices,
        #        "indptr":tmp.indptr,
        #        "geneIds":combined["geneIds"],
        #        "sampleIds":combined["sampleIds"]
        #        }

        #    output["chunk" + str(i+1)] = assay
        #    del tmp,assay
        #    gc.collect()

        # Store back as Dense Form
        print("FINISH imputing one chunk!",flush=True)
        new_assay = {
                "matrix":imputed.T.to_numpy().tolist(),
                "geneIds": combined["geneIds"],
                "sampleIds": combined["sampleIds"]
                }
        output["chunk"+str(i+1)] = compress_assay(new_assay)
        print("FINISH storing one chunk!", flush=True)


        nb_genes += len(set(model.targets.flatten()))
        r, _ = model.score(frameddata)
        sum_r += r 
        rows, cols = frameddata.shape
        data_dropout += calc_dropout(data) * 100
        impu_dropout += calc_dropout(imputed.to_numpy()) * 100

        if i == len(chunks)-1:
            fig, ax = plt.subplots(1, 2)
            LABELS_PARAMS = {"fontsize": 14, "fontweight": "bold", "fontname": "serif"}
            
            
            vmax = np.percentile(np.log10(1 + data.flatten()), 99)
            print("Generating Heatmap", flush=True)
            ax[0].imshow(np.log10(1 + data), aspect="auto", vmax=vmax)
            ax[1].imshow(np.log10(1 + imputed), aspect="auto", vmax=vmax)
            ax[0].set_xlabel("Genes", **LABELS_PARAMS)
            ax[1].set_xlabel("Genes", **LABELS_PARAMS)
            ax[0].set_ylabel("Cells", **LABELS_PARAMS)
            ax[0].set_title("last chunk raw (log)", **LABELS_PARAMS)
            ax[1].set_title("last chunk imputed (log)", **LABELS_PARAMS)
            gn.add_current_figure_to_results("Heatmaps")

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
                output["origin data size"][1],
                output["origin data size"][0],
                nb_genes,
                data_dropout/(len(chunks)),
                impu_dropout/(len(chunks)),
                sum_r/(len(chunks))
            )
            gn.add_result(message, data_type="markdown")
            del combined,data,frameddata,imputed, new_assay
            gc.collect()
            break
        
        del combined,data,frameddata,imputed, new_assay
        gc.collect()

    #for i in range(len(output)-1):
    #    if i == 0:
    #        output["chunk"+str(i+1)]["data"] =  output["chunk"+str(i+1)]["data"].tolist()
    #        output["chunk"+str(i+1)]["indices"] =  output["chunk"+str(i+1)]["indices"].tolist()
    #        output["chunk"+str(i+1)]["indptr"] =  output["chunk"+str(i+1)]["indptr"].tolist()

    gn.export_statically(output, "Imputed assay")
    gn.commit()

    #print("--- %s seconds --- for compress" % (time.time() - deepimpute_time), flush = True)
    #print("--- %s seconds --- for whole gbox" % (time.time() - start_time), flush = True)
    #time.sleep(10)

if __name__ == "__main__":
    main()
