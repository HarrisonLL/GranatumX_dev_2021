import math

import scanpy.api as sc
import numpy as np
from granatum_sdk import Granatum
from granatum_sdk_subclass import granatum_extended
import time

def main():
    start_time = time.time()
    print("Start to switch", flush = True)
    gn = granatum_extended("gbox-scanpygenefilering")
    output = {}

    chunks = gn.get_import('assay')
    new_chunks = gn.adjust_transform_test(chunks)
    switch_time = time.time()
    print("--- %s seconds --- for transforming" % (switch_time - start_time), flush = True)

    chunk1_assay = gn.combine_new_chunk(new_chunks["chunk1"])

    adata = gn.ann_data_from_assay(chunk1_assay)
    min_cells_expressed = gn.get_arg("min_cells_expressed")
    min_mean = gn.get_arg("min_mean")
    max_mean = gn.get_arg("max_mean")
    min_disp = gn.get_arg("min_disp")
    max_disp = gn.get_arg("max_disp")

    num_genes_before = adata.shape[1]
    
    print("Start to filter", flush = True)

    sc.pp.filter_genes(adata, min_cells=min_cells_expressed)

    filter_result = sc.pp.filter_genes_dispersion(
        adata.X, flavor='seurat', min_mean=math.log(min_mean), max_mean=math.log(max_mean), min_disp=min_disp, max_disp=max_disp,
    )
    adata = adata[:, filter_result.gene_subset]

    sc.pl.filter_genes_dispersion(filter_result)

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
                "Number of genes after filtering: **{}**".format(adata.shape[1]),
            ]
        ),
        type="markdown",
    )
    
    print("Compressing", flush = True)
    output["chunk1"] = gn.compress_chunk(gn.assay_from_ann_data(adata))
    gn.export(output, "Filtered Assay", dynamic=False)

    gn.commit()
    print("--- %s seconds --- for whole gbox" % (time.time() - start_time), flush = True)
    time.sleep(20)

if __name__ == "__main__":
    main()
