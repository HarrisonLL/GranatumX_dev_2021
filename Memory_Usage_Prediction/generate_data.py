import numpy as np
import pandas as pd
import os
import shutil

def generate_data_matrix(cell_sizes, gene_sizes, percent_nz):
    datasets = {}
    for percent in percent_nz:
        percent = percent / 100
        for csize in cell_sizes:
            for gsize in gene_sizes:
                matrix = np.zeros((csize*gsize,))
                indices = np.random.choice(csize*gsize,size=int(csize*gsize*percent),replace=False)
                for i in indices:
                    matrix[i] = np.random.randint(low=1,high=10000) # FIX ME
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
    gene_sizes = [100,500,1000,3000,10000]
    percent_nz = [1,5,10,20] 
    datasets = generate_data_matrix(cell_sizes,gene_sizes,percent_nz)
    for k,v in datasets.items():
        csv_name = "_".join(k)+".csv"
        pd.DataFrame(data=v).to_csv(os.path.join(dirpath,csv_name))




if __name__ == "__main__":
    main()
