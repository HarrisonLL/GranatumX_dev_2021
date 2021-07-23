import numpy as np
import pandas as pd
import os
import shutil
from tqdm import tqdm
from scipy import sparse
from scipy.sparse import rand
import gc


def main():
     # destroy if exits, and then create ./datasets directory
    dirpath = './datasets'
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
    os.mkdir(dirpath)

    cell_sizes = [100, 500,1000, 2000, 4000,6000,8000]
    gene_sizes = [int(1e4), int(2e4), int(3e4), int(4e4), int(5e4),int(6e4)]
    percent_nz = [0.1,0.5,1,2,5,10,20,30]
    print("Generating datasets and saving them..", flush=True)
    for csize in tqdm(cell_sizes):
        for gsize in tqdm(gene_sizes,leave=False):
            for percent in tqdm(percent_nz, leave=False):
                percent /= 100
                matrix = rand(gsize, csize, density=percent, dtype=np.uint8, format="csr", random_state=42)
                name = "_".join((str(gsize), str(csize), str(percent*100))) + ".npz"
                sparse.save_npz(os.path.join('./datasets',name), matrix)
                del (matrix)
                gc.collect()

if __name__ == "__main__":
    main()
