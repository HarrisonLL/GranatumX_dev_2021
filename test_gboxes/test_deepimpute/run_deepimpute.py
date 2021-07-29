import granatum_sdk
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from deepimpute.deepImpute import deepImpute 
import time
import json
import gzip
from io import StringIO, BytesIO
import base64
import gc

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
    gn = granatum_sdk.Granatum()
    output = {}
    
    #firstly save to chunks
    chunks = gn.get_import("assay")
    print(len(chunks.keys()), flush = True)

    #assay = gn.get_import("assay")
    assay = decode_json(chunks["chunk1"])
    data = np.array(assay.get("matrix")).T
    barcodes = assay.get("sampleIds")
    genes = assay.get("geneIds")
    del assay
    gc.collect()

    seed = gn.get_arg("seed")
    checkbox = gn.get_arg("use_auto_limit")
    cell_subset = gn.get_arg("cell_subset")

    NN_lim = {False: gn.get_arg("NN_lim"), True: "auto"}.get(checkbox, True)

    #model = MultiNet(n_cores="all", seed=seed)
    #model.fit(data, NN_lim=NN_lim, cell_subset=cell_subset)

    frameddata = pd.DataFrame(data)
    print("imputed", flush = True)
    imputed, model = deepImpute(frameddata, NN_lim=NN_lim, cell_subset=cell_subset)

    LABELS_PARAMS = {"fontsize": 14, "fontweight": "bold", "fontname": "serif"}
    print("Calculate vmax", flush = True)
    
    #tmp data, imputed
    data = data + 1
    imputed = imputed + 1

    vmax = np.percentile(np.log10(data.flatten()), 99)

    print("Generating Heatmap", flush = True)
    fig, ax = plt.subplots(1, 2)
    ax[0].imshow(np.log10(data), aspect="auto", vmax=vmax)
    ax[1].imshow(np.log10(imputed), aspect="auto", vmax=vmax)
    ax[0].set_xlabel("Genes", **LABELS_PARAMS)
    ax[1].set_xlabel("Genes", **LABELS_PARAMS)
    ax[0].set_ylabel("Cells", **LABELS_PARAMS)
    ax[0].set_title("raw (log)", **LABELS_PARAMS)
    ax[1].set_title("imputed (log)", **LABELS_PARAMS)
    
    print("Heatmap created!", flush = True)

    gn.add_current_figure_to_results("Heatmaps")
    nb_genes = len(set(model.targets.flatten()))
    #nb_genes = np.sum([len(net.targetGenes) for net in model.targets])
    print("finished adding", flush = True)
    del data
    gc.collect()
    def calc_dropout(matrix):
        return np.sum((np.array(matrix) == 0)) * 1. / data.size

   # r, p = model.score(frameddata)
   # rows, cols = frameddata.shape
    print("start to transform", flush = True)
    #transform data, imputed back
    #data = data - 1
    imputed = imputed - 1
    #print("creating message", flush = True)

    #message = "\n".join(
    #    [
    #        "  - Data frame number of rows: **{0}**",
    #        "  - Data frame number of columns: **{1}**",
    #        "  - Number of imputed genes: **{2}**",
    #        "  - Percentage of dropout entries *before* imputation: **{3:.2f}%**",
    #        "  - Percentage of dropout entries *after* imputation: **{4:.2f}%**",
    #        "  - Accuracy (correlation) on masked data: **{5:.2f}**"
    #    ]
    #).format(
    #    rows,
    #    cols,
    #    nb_genes,
    #    calc_dropout(data) * 100,
    #    calc_dropout(imputed.to_numpy()) * 100,
    #   r
    #)

    #gn.add_result(message, data_type="markdown")

    exported_assay = {
                    "matrix":  imputed.T.to_numpy().tolist(),
                    "sampleIds": barcodes,
                    "geneIds": genes
                }
    output["chunk1"] = compress_assay(exported_assay)
    #gn.export_statically(assay, "Imputed assay")
    with open(os.path.join(gn.exports_dir,"chunks"),"wt") as f:
        json.dump(output, f)
    gn.dynamic_exports.append({"extractFrom": "chunks", "kind": "assay", "meta": None})
    gn.commit()
    del imputed, exported_assay
    gc.collect()
    print("--- %s seconds ---" % (time.time() - start_time), flush = True)
    time.sleep(10)

if __name__ == "__main__":
    main()
