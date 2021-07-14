from deepimpute.multinet import MultiNet
import pandas as pd
import numpy as np
import psutil
import os
import tracemalloc
import csv
from tqdm import tqdm
np.random.seed(12345)


def get_ava_memory():
    usage = psutil.virtual_memory()
    ava_mb = usage.available/(1024*1024)
    print("Ava memory (MB):", ava_mb, flush=True)
    return ava_mb


def run_deep_impute(data):
    # Using custom parameters
    print("Start to Fit and Predict ==============>", flush=True)
    NN_params = {
            'learning_rate': 1e-4,
            'batch_size': 64,
            'max_epochs': 200,
            'ncores': 5,
            'sub_outputdim': 512,
            'architecture': [
                {"type": "dense", "activation": "relu", "neurons": 200},
                {"type": "dropout", "activation": "dropout", "rate": 0.3}]
        }
    multinet = MultiNet(**NN_params, verbose=None)
    multinet.fit(data,cell_subset=1,minVMR=0.5)
    imputedData = multinet.predict(data)
    print("============> End Fit and Predict", flush=True)


def main():
    ava_mb = get_ava_memory()
    memory_data = []
    # read data from ./datasets and perform deepimpute
    # save the result to rows
    rows = []
    for file in tqdm(os.listdir("./datasets")):
        data = pd.read_csv(os.path.join("./datasets", file))
        data.drop(columns=["Unnamed: 0"],inplace=True)
        tracemalloc.start()
        run_deep_impute(data.T)
        _,PEAK = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        memory_data.append(PEAK/1024/1024)
        print("Peak Usage of the function:" + str(PEAK/1024/1024),flush=True)
        tmp = file.replace(".csv", "").split("_")
        gene_size = tmp[0]
        cell_size = tmp[1]
        percent = tmp[2]
        rows.append([int(gene_size), int(cell_size),float( percent), PEAK])
        del(data)
        del(PEAK)
        del(tmp)
        del(gene_size)
        del(cell_size)
        del(percent)
        
    # write the data to csv
    performance = pd.DataFrame(data=rows, columns=["Gene Size", "Cell Size", "Percent", "Peak Memory Usage"])
    performance.to_csv("performance.csv")


if __name__ == "__main__":
    main()
