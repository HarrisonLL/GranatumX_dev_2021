#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time

from granatum_sdk import Granatum
from granatum_sdk_subclass import granatum_extended

def main():
    start_time = time.time()
    gn = granatum_extended("gbox-logtransformation")
    output = {}

    chunks = gn.get_import('assay')
    assay = gn.decompress_chunk(chunks["chunk1"])
    #assay = gn.get_import('assay')
    matrix = np.array(assay.get('matrix'))

    log_base = gn.get_arg('logBase')
    pseudo_counts = gn.get_arg('pseudoCounts')

    transformed_matrix = np.log(matrix + pseudo_counts) / np.log(log_base)
    print("point 0", flush = True)
    non_zero_values_before = matrix.flatten()
    # non_zero_values_before = non_zero_values_before[(non_zero_values_before > np.percentile(non_zero_values_before, 5)) &
    #                                                 (non_zero_values_before < np.percentile(non_zero_values_before, 95))]
    print("point 1", flush = True)
    non_zero_values_before = non_zero_values_before[(non_zero_values_before > np.percentile(non_zero_values_before, 5))]
    print("point 2", flush = True)
    non_zero_values_after = transformed_matrix.flatten()
    # non_zero_values_after = non_zero_values_after[(non_zero_values_after > np.percentile(non_zero_values_after, 5)) &
    #                                               (non_zero_values_after < np.percentile(non_zero_values_after, 95))]
    print("point 3", flush = True)
    non_zero_values_after = non_zero_values_after[(non_zero_values_after > np.percentile(non_zero_values_after, 5))]

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

    assay['matrix'] = transformed_matrix.tolist()
    output["chunk1"] = gn.compress_chunk(assay)
    gn.export_statically(output, 'Log transformed assay')

    gn.commit()
    print("--- %s seconds ---" % (time.time() - start_time), flush = True)
    time.sleep(10)


if __name__ == '__main__':
    main()
