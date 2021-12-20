#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time
import scipy

from granatum_sdk import Granatum
from granatum_sdk_subclass import granatum_extended
from granatum_sdk_subclass2 import granatum_extended2

def main():
    start_time = time.time()
    bool_sparse = False
    gn = granatum_extended("logtransformation")

    chunks = gn.get_import('assay')
    print(chunks["current chunk"], flush = True)
    if chunks["current chunk"][-1] == "sparse":
        print("it's sparse", flush = True)
        gn = granatum_extended2("logtransformation")
        chunks = gn.get_import('assay')
        bool_sparse = True

    output = {"origin data size":chunks["origin data size"], 
              "current chunk":["logtransformation", "col", "sparse"],
              "suggested chunk":chunks["suggested chunk"]}

    chunks["suggested chunk"]["logtransformation"] = ["col", 3456] 
    print("Start to switch", flush = True)
    gn.adjust_transform(chunks)
    switch_time = time.time()
    print("--- %s seconds --- for transforming" % (switch_time - start_time), flush = True)
    chunks = gn.new_file

    log_base = gn.get_arg('logBase')
    pseudo_counts = gn.get_arg('pseudoCounts')

    for i in range(len(chunks)):
        print("current working on chunk %s" % (i + 1), flush = True)
        combined = gn.combine_new_chunk(chunks["chunk"+str(i+1)], "col")
        if bool_sparse:
            matrix  = scipy.sparse.csc_matrix((combined.get("data"), combined.get("indices"), combined.get("indptr")), shape=(len(combined["geneIds"]), len(combined["sampleIds"]))).todense()
            combined#["matrix"] = matrix
            del(combined["data"])
            del(combined["indices"])
            del(combined["indptr"])
    #assay = gn.decompress_chunk(chunks["chunk1"])
    #assay = gn.get_import('assay')
    #matrix = np.array(assay.get('matrix'))
        else:
            matrix = combined.get('matrix')

        transformed_matrix = np.log(matrix + pseudo_counts) / np.log(log_base)
    #print("point 0", flush = True)
        non_zero_values_before = matrix.flatten()
    # non_zero_values_before = non_zero_values_before[(non_zero_values_before > np.percentile(non_zero_values_before, 5)) &
    #                                                 (non_zero_values_before < np.percentile(non_zero_values_before, 95))]
    #print("point 1", flush = True)
        non_zero_values_before = non_zero_values_before[(non_zero_values_before > np.percentile(non_zero_values_before, 5))]
    #print("point 2", flush = True)
        non_zero_values_after = transformed_matrix.flatten()
    # non_zero_values_after = non_zero_values_after[(non_zero_values_after > np.percentile(non_zero_values_after, 5)) &
    #                                               (non_zero_values_after < np.percentile(non_zero_values_after, 95))]
    #print("point 3", flush = True)
        non_zero_values_after = non_zero_values_after[(non_zero_values_after > np.percentile(non_zero_values_after, 5))]

        if i == 0:
            plt.figure()

            plt.subplot(2, 1, 1)
            plt.title('Before log transformation')
            plt.hist(non_zero_values_before, bins=100)
            plt.ylabel('Frequency')
            plt.xlabel('Expression level')

            plt.subplot(2, 1, 2)
            plt.title('After log transformation')
            plt.hist(non_zero_values_after, bins=100)
            plt.ylabel('Frequency')
            plt.xlabel('Expression level')

            plt.tight_layout()
            print("Finished plotting", flush = True)
            caption = (
                'The distribution of expression level before and after log transformation. Only the values greater '
                'than the 5 percentile (usually zero in single-cell data) and lower than 95 percentile are considered.'
            )
            gn.add_current_figure_to_results(caption, zoom=2, dpi=50)

        combined['matrix'] = transformed_matrix.tolist()
        output["chunk"+str(i+1)] = gn.compress_chunk(combined)
    gn.export(output, 'Log transformed assay', dynamic = False)
    gn.commit()
    print("--- %s seconds ---" % (time.time() - start_time), flush = True)
    time.sleep(10)


if __name__ == '__main__':
    main()
