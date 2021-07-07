from deepimpute.multinet import MultiNet
from sklearn.linear_model import LinearRegression
import loompy
import pandas as pd
import numpy as np
import argparse
import psutil
from resource import getrusage, RUSAGE_SELF
import matplotlib.pyplot as plt


def get_ava_memory():
    usage = psutil.virtual_memory()
    ava_mb = usage[1]/(1024*1024)
    print("Ava memory (MiB):", ava_mb, flush=True)
    return ava_mb


def get_peak_memory():
    PEAK =  int(getrusage(RUSAGE_SELF).ru_maxrss / (1024*1024))
    print("Peak memory (MiB):", PEAK, flush=True)
    return PEAK


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
    PEAK = get_peak_memory()
    return PEAK


def main(argv):
    file_path = argv[0]
    gene_size = argv[1]
    data = loompy.connect(file_path)
    ava_mb = get_ava_memory()
    memory_data = []
    cell_sizes = [500, 1000, 5000, 10000, 15000, 20000]
    for size in cell_sizes:
        df = pd.DataFrame((data[:gene_size, :size].T))
        PEAK = run_deep_impute(df)
        memory_data.append(PEAK)
    # run linear regression to get coeff
    reg = LinearRegression().fit(np.array(cell_sizes).reshape(6,1), np.array(memory_data).reshape(6,))
    coef = reg.coef_[0]
    inter = reg.intercept_
    print("Suggested cell size:", (ava_mb-inter)//coef)
    plt.figure(figsize=(4, 3))
    ax = plt.axes()
    ax.scatter(cell_sizes, memory_data)
    ax.plot(np.linspace(0, 22000, 220000), reg.predict(np.linspace(0, 22000, 220000).reshape(-1, 1)))
    plt.savefig("linear_regression.png")





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file_path')
    parser.add_argument('--gene_size')
    args = parser.parse_args()
    argv = [None]*2
    for arg in vars(args):
        value = getattr(args, arg)
        if arg == "input_file_path":
            argv[0] = value
        elif arg == "gene_size":
            argv[1] = int(value)
    main(argv)