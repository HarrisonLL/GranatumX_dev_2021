import numpy as np
import pandas as pd
import os
import shutil
from tqdm import tqdm

def generate_data_matrix(cell_sizes, gene_sizes, percent_nz):
    datasets = {}
    for csize in tqdm(cell_sizes):
        for gsize in gene_sizes:
            for percent in percent_nz:
                percent = percent / 100
                matrix = np.zeros((csize*gsize,))
                indices = np.random.choice(csize*gsize,size=int(csize*gsize*percent),replace=False)
                for i in indices:
                    matrix[i] = np.random.randint(low=1,high=100) # FIX ME
                matrix = matrix.reshape((gsize, csize))
                datasets[(str(gsize), str(csize), str(percent*100))] = matrix
    return datasets

def main():
     # destroy if exits, and then create ./datasets directory
    dirpath = './datasets'
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
    os.mkdir(dirpath)

    cell_sizes = [1000,2000,5000,10000,20000]
    gene_sizes = [100,500,1000,5000]
    percent_nz = [5,10,20,40]
    print("Generating datasets..", flush=True)
    datasets = generate_data_matrix(cell_sizes,gene_sizes,percent_nz)
    print("Saving..",flush=True)
    for k,v in tqdm(datasets.items()):
        csv_name = "_".join(k)+".csv"
        pd.DataFrame(data=v).to_csv(os.path.join(dirpath,csv_name))




if __name__ == "__main__":
    main()
