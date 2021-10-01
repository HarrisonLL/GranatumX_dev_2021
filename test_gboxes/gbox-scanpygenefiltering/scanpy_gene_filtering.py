import math
import gc
import scipy
import scanpy.api as sc
import numpy as np
from granatum_sdk_subclass import granatum_extended
from granatum_sdk_subclass2 import granatum_extended2
import time

def main():
    
    start_time = time.time()
    bool_sparse = False
    gn = granatum_extended("scanpygenefilering")
    chunks = gn.get_import('assay')
    org_chunk_kind = chunks["current chunk"][-2]

    if chunks["current chunk"][-1] == "sparse":
        gn = granatum_extended2("scanpygenefilering")
        chunks = gn.get_import('assay')
        bool_sparse = True
    output = {"origin data size":chunks["origin data size"], "current chunk":["scanpygenefilering", "row", "dense"], "suggested chunk":chunks["suggested chunk"]}
    
    print("Start to switch", flush = True)
    gn.adjust_transform(chunks)
    switch_time = time.time()
    print("--- %s seconds --- for transforming" % (switch_time - start_time), flush = True)
    chunks = gn.new_file

    num_genes_before = 0
    num_genes_after = 0

    for i in range(len(chunks)):
        combined = gn.combine_new_chunk(chunks["chunk"+str(i+1)], org_chunk_kind)
        if bool_sparse:
            matrix  = scipy.sparse.csc_matrix((combined.get("data"), combined.get("indices"), combined.get("indptr")), shape=(len(combined["geneIds"]), len(combined["sampleIds"]))).todense()
            combined["matrix"] = matrix
            del(combined["data"])
            del(combined["indices"])
            del(combined["indptr"])

        adata = gn.ann_data_from_assay(combined)
        min_cells_expressed = gn.get_arg("min_cells_expressed")
        min_mean = gn.get_arg("min_mean")
        max_mean = gn.get_arg("max_mean")
        min_disp = gn.get_arg("min_disp")
        max_disp = gn.get_arg("max_disp")

        num_genes_before += adata.shape[1]
        
        print("Start to filter", flush = True)

        sc.pp.filter_genes(adata, min_cells=min_cells_expressed)

        filter_result = sc.pp.filter_genes_dispersion(
            adata.X, flavor='seurat', min_mean=math.log(min_mean), max_mean=math.log(max_mean), min_disp=min_disp, max_disp=max_disp,
        )
        adata = adata[:, filter_result.gene_subset]

        sc.pl.filter_genes_dispersion(filter_result)
        
        num_genes_after += adata.shape[1]

        
        print("Compressing", flush = True)
        output["chunk"+str(i+1)] = gn.compress_chunk(gn.assay_from_ann_data(adata))
    

        if i == len(chunks) - 1:
            print("Add to result", flush = True)
            gn.add_current_figure_to_results(
                "Each dot represent a gene. The gray dots are the removed genes. The x-axis is log-transformed.",
                zoom=3,
                dpi=50,
                height=400,
            )

            gn.add_result(
                "\n".join(
                [
                "Number of genes before filtering: **{}**".format(num_genes_before),
                "",
                "Number of genes after filtering: **{}**".format(num_genes_after),
                ]
                ),
             type="markdown",
            )
        del adata, min_cells_expressed, min_mean, max_mean, min_disp, max_disp,filter_result
        gc.collect()
    
    gn.export(output, "Filtered Assay", dynamic=False)
    gn.commit()
    print("--- %s seconds --- for whole gbox" % (time.time() - start_time), flush = True)
    time.sleep(20)

if __name__ == "__main__":
    main()
