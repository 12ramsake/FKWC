# FKWC
Multisample testing for covariance operator difference in R

See 'Robust nonparametric hypothesis tests for differences in the covariance structure of functional data':  https://arxiv.org/abs/2106.10173


Note that the code of other authors to generate p-values from other test statistics is missing, as it is not mine to distribute.


I also do not include the financial data analysis, since the data is licensed


If you just want to use an implementation of the test, just choose a package that computes the depth values you want to use, e.g., fda.usc, and compute the depth values
After you have the depth values, just use the rank test function in R; kruskal.wallis(). If you have trouble running the test, please contact me at kramsay2@yorku.ca
