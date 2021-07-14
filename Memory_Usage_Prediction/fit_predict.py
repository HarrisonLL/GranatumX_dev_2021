import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import cross_val_score


"""fit linear regression, quadratic linear regression, cubic linear regression
return cross validaton r square score for each model"""
def model_fitting(X,y):
    reg = LinearRegression()
    rs1 = np.mean(cross_val_score(reg, X, y, cv=5))

    poly = PolynomialFeatures(degree=2)
    X_new2 = poly.fit_transform(X)
    reg = LinearRegression()
    rs2 = np.mean(cross_val_score(reg, X_new2, y, cv=5))

    poly = PolynomialFeatures(degree=3)
    X_new3 = poly.fit_transform(X)
    reg = LinearRegression()
    rs3 = np.mean(cross_val_score(reg, X_new3, y, cv=5))
   
    X_new9 = [[1,x[0],x[0]**2, x[1],x[1]**2, x[2], x[2]**2] for x in X]
    X_new10 = [[1, 2*x[0]*x[1], 2*x[1]*x[2]] for x in X]

    rs = [rs1, rs2, rs3]
    for X_new in [ X_new9, X_new10]:
        reg = LinearRegression()
        temp = np.mean(cross_val_score(reg, X_new, y, cv=5))
        rs.append(temp)

    return rs

def main():
    df = pd.read_csv("performance.csv")
    rs = model_fitting(df.drop(columns=["Peak Memory Usage"]).to_numpy(), df["Peak Memory Usage"].to_numpy())
    print(rs)
    print(np.argmax(rs))


if __name__ == "__main__":
    main()

