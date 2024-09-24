# Repo `yieldcurve`
It is mainly focus on yield curve in financial markets.

Yield curves are important, especially in derivatives world. Sometime, we are going to convert dirty prices into yields and then fit discrete yields into a continuous curve. This repositoey contains project code that can compute yields from prices and methods to create a yield curve, using Statisticality. May be updated in the future.

# Current methods
* version 240923
1) Linear Regression using OLS, with formula like `y = a + bx + cx^2 + dx^3 + f*exp(-x) + g*log(x) + h*x*exp(-x)`
2) Non-linear Regression using Gradient Descending, with formula like `y = a + bx + cx^2 + dx^3 + f*exp(g*x) + h*log(x+1) + j*x*exp(k*x)`
3) Neural Network from `draft` with feature engineering (non-parametic, mostly precise)
4) utils: `drop_outliers`, `bond_pricing`, ...
