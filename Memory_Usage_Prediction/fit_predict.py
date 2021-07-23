import pandas as pd
import numpy as np
import sklearn
from sklearn.linear_model import LinearRegression
from sklearn import linear_model
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import cross_val_score
from sklearn.metrics import mean_squared_error
from sympy import symbols, Eq, solve
import psutil

"""fit linear regression, quadratic linear regression, cubic linear regression
return cross validaton r square, mse for each model, coeff for the best model"""
def model_fitting(X,y):
    poly = PolynomialFeatures(degree=2)
    X = poly.fit_transform(X)
    reg = linear_model.Lasso(alpha=1.0,max_iter=50, selection="random",warm_start=False,random_state=0)
    rs = np.mean(cross_val_score(reg,X,y,cv=10, scoring="r2"))
    rmse = np.mean(np.sqrt(np.absolute(cross_val_score(reg,X,y,scoring='neg_mean_squared_error',cv=10))))
    return rs, rmse, reg.fit(X,y).coef_


def pred_cell_size(coef,genesize, percent,ava_mem):
    x = symbols('x')
    eq1 = Eq(coef[0] + coef[1]*genesize + coef[2]*x + coef[3]*percent + coef[4]*(genesize**2) + coef[5]*(genesize*x) + coef[6]*(genesize*percent)+ coef[7]*(x**2) + coef[8]*(x*percent)+coef[9]*(percent**2)-ava_mem)
    return [int(sol) for sol in solve(eq1) if sol>0]


def model_time(X,y):
    poly = PolynomialFeatures(degree=2)
    X_new2 = poly.fit_transform(X)
    #reg = LinearRegression(fit_intercept=False)
    reg = linear_model.Lasso(alpha=1.0,max_iter=50, selection="random",warm_start=False,random_state=0)
    rs2 = np.mean(cross_val_score(reg, X_new2, y, cv=10))
    mse2 = np.mean(np.sqrt(np.absolute(cross_val_score(reg,X,y,scoring='neg_mean_squared_error',cv=10))))
    return rs2, mse2, reg.fit(X_new2,y).coef_


def predict_time(coef,genesize,cellsize,percent):
    sol = coef[0] + coef[1]*genesize + coef[2]*cellsize + coef[3]*percent + coef[4]*(genesize**2) + coef[5]*(genesize*cellsize) + coef[6]*(genesize*percent)+ coef[7]*(cellsize**2) + coef[8]*(cellsize*percent)+coef[9]*(percent**2)
    return sol


def pred_cell_size(coef,genesize, percent,ava_mem):
    x = symbols('x')
    eq1 = Eq(coef[0] + coef[1]*genesize + coef[2]*x + coef[3]*percent + coef[4]*(genesize**2) + coef[5]*(genesize*x) + coef[6]*(genesize*percent)+ coef[7]*(x**2) + coef[8]*(x*percent)+coef[9]*(percent**2)-ava_mem)
    return [int(sol) for sol in solve(eq1) if sol>0]



def main():
    df = pd.read_csv("performance.csv")
    rs, mse, coef = model_fitting(df[["Gene Size", "Cell Size","Percent"]].to_numpy(), df["Peak Memory Usage"].to_numpy())
    print(rs, mse, coef)
    df2 = pd.read_csv("time_performance.csv")
    rs2, mse2, coef2 = model_time(df2[["Gene Size", "Cell Size","Percent"]].to_numpy(), df2["Time"].to_numpy())
    print(rs2,mse2,coef2)
    #ava_mem = (psutil.virtual_memory()).available/(1024*1024)
    #print("Current available memory is %.3f GB"%(ava_mem/1024))
    #print("If the gene size is 20000, and percent of non-zeros is 2%, we can chunk cell size to")
    #print()
    #cellsize = pred_cell_size(coef, 20000, 2, ava_mem)[0]
    #print(cellsize)
    
    #df2 = pd.read_csv("time_performance.csv")
    #df2 = df2.drop(columns=["Unnamed: 0"])
    #rs, mse, coef2 = model_time(df2.drop(columns=["Time"]).to_numpy(), df2["Time"].to_numpy())
    #print("It will takes %.3f minutes" %(predict_time(coef2,8000,cellsize,10)/60)) 
    pd.DataFrame(data=[coef,coef2], columns = ["intercept", "gene", "cell", "percent", "gene**2", "gene*cell", "gene*percent", "cell**2", "cell*percet", "percent**2"], index = ["memory", "time"]).to_csv("coeffs.csv")
    #print(pred_cell_size(coef, 8000, 10, ava_mem))
    #print()
    #print(mse2, mse)

if __name__ == "__main__":
    main()

