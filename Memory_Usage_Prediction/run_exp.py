from deepimpute.multinet import MultiNet
import pandas as pd
import numpy as np
from scipy import sparse
import psutil
import os
import tracemalloc
import csv
from tqdm import tqdm
from time import time
import gc



def get_ava_memory():
    usage = psutil.virtual_memory()
    ava_mb = usage.available/(1024*1024)
    print("Ava memory (MB):", ava_mb, flush=True)



def run_deep_impute(data):
    # Simulate DeepImpute Gbox
    print("Start to Fit and Predict ==============>", flush=True)
    np.random.seed(12345)
    data = pd.DataFrame(data).T
    multi = MultiNet()
    multi.fit(data, NN_lim="auto", cell_subset=1, minVMR=0.5)
    multi.predict(data, imputed_only=False, policy="restore")
    print("============> End Fit and Predict", flush=True)


def main():

    # read data from ./datasets and perform deepimpute
    # save the result to rows
    with open("performance3.csv", "w") as f:
        with open("time_performance3.csv", "w") as f2:
            f.write("Gene Size,Cell Size,Percent,Peak Memory Usage\n")
            f2.write("Gene Size,Cell Size,Percent,Time\n")
            count = 0
            for file in tqdm(sorted(os.listdir("./datasets"))):
                if 152 <= count < 200:
                    print(file,flush=True)
                    get_ava_memory()
                    file_path = os.path.join("./datasets", file)
                    data = sparse.load_npz(file_path).todense()
                    print(data.shape,flush=True)
                
                    begin = time()
                    tracemalloc.start()
                    run_deep_impute(data)
                    _,PEAK = tracemalloc.get_traced_memory()
                    tracemalloc.stop()
                    end = time()
                    elapse = end - begin

                    print("Peak Usage of the function:" + str(PEAK/1024/1024),flush=True)
                    tmp = file.replace(".npz", "").split("_")
                    gene_size = tmp[0]
                    cell_size = tmp[1]
                    percent = tmp[2]
                    f.write(",".join([gene_size, cell_size, percent, str(PEAK/1024/1024)])+ "\n")
                    f2.write(",".join([gene_size, cell_size, percent, str(elapse)])+ "\n")
                    del(data)
                    del(PEAK)
                    del(tmp)
                    del(gene_size)
                    del(cell_size)
                    del(percent)
                    del(begin)
                    del(end)
                    del(elapse)
                    gc.collect()
                count += 1
    f.close()
    f2.close()


if __name__ == "__main__":
    main()
