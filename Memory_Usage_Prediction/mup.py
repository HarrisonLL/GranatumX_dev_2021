from deepimpute.multinet import MultiNet
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import mean_squared_error
import loompy
import pandas as pd
import numpy as np
import argparse
import psutil
from resource import getrusage, RUSAGE_SELF
import matplotlib.pyplot as plt
np.random.seed(12345)


def get_ava_memory():
    usage = psutil.virtual_memory()
    ava_mb = usage.available/(1024*1024)
    print("Ava memory (MB):", ava_mb, flush=True)
    return ava_mb


def get_peak_memory():
    PEAK =  int(getrusage(RUSAGE_SELF).ru_maxrss * 1.04858 / (1024)) # MiB to MB
    print("Peak memory (MB):", PEAK, flush=True)
    return PEAK


def generate_data_matrix(cell_sizes, gene_sizes, percent_nz):
    datasets = {}
    for percent in percent_nz:
        percent = percent / 100
        for csize in cell_sizes:
            for gsize in gene_sizes:
                matrix = np.zeros((csize*gsize,))
                indices = np.random.randint(0,csize*gsize,size=int(csize*gsize*percent))
                for i in indices:
                    matrix[i] = np.random.randint(0,10) # FIX ME
                matrix = matrix.reshape((gsize, csize))
                datasets[(gsize, csize, percent*100)] = matrix
    return datasets
    

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

"""fit linear regression, quadratic linear regression, cubic linear regression
return mse for each model"""
def model_fitting(X,y):
    X_new1 = [[x[0], x[1], x[2]] for x in X]
    reg = LinearRegression().fit(X_new1, y)
    y_pred1 = reg.predict(X_new1)
    mse1 = mean_squared_error(y, y_pred1)

    poly = PolynomialFeatures(degree=2)
    X_new2 = poly.fit_transform(X_new1)
    y_pred2 = reg.predict(X_new2)
    mse2 = mean_squared_error(y, y_pred2)

    poly = PolynomialFeatures(degree=3)
    X_new3 = poly.fit_transform(X_new1)
    y_pred3 = reg.predict(X_new3)
    mse3 = mean_squared_error(y, y_pred3)
    return mse1, mse2, mse3


def main():
    
    ava_mb = get_ava_memory()
    memory_data = []
    cell_sizes = [500, 1000, 5000, 10000, 15000, 20000]
    gene_sizes = [500, 1000, 5000, 10000, 15000, 20000]
    percent_nz = [0.01, 0.1, 0.5, 1, 10, 20, 50]
    datasets = generate_data_matrix(cell_sizes, gene_sizes, percent_nz)
    for data in datasets.values():
        PEAK = run_deep_impute(pd.DataFrame(data.T))
        memory_data.append(PEAK)
    mse1, mse2, mse3 = model_fitting(datasets.keys(), memory_data)
    print(mse1, mse2, mse3)

    # run linear regression to get coeff
    # reg = LinearRegression().fit(np.array(cell_sizes).reshape(6,1), np.array(memory_data).reshape(6,))
    # coef = reg.coef_[0]
    # inter = reg.intercept_
    # print("Suggested cell size:", (ava_mb-inter)//coef)
    # plt.figure(figsize=(4, 3))
    # ax = plt.axes()
    # ax.scatter(cell_sizes, memory_data)
    # ax.plot(np.linspace(0, 22000, 220000), reg.predict(np.linspace(0, 22000, 220000).reshape(-1, 1)))
    # plt.savefig("linear_regression.png")





if __name__ == "__main__":
    main()