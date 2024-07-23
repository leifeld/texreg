# texreg

Conversion of R Regression Output to LaTeX or HTML Tables.

Converts coefficients, standard errors, significance stars, and goodness-of-fit statistics of statistical models into LaTeX tables or HTML tables/MS Word documents or to nicely formatted screen output for the R console for easy model comparison. A list of several models can be combined in a single table. The output is highly customizable. New model types can be easily implemented.

## Documentation

Details on **texreg** can be found in the following article:

Leifeld, Philip (2013): **texreg**: Conversion of Statistical Model Output in R to LaTeX and HTML Tables. Journal of Statistical Software 55(8): 1-24. doi:[10.18637/jss.v055.i08](http://dx.doi.org/10.18637/jss.v055.i08)

An updated version of this paper is included as a [package vignette](https://cran.r-project.org/web/packages/texreg/vignettes/texreg.pdf). Additional details on how to write custom extensions for **texreg** and manipulate **texreg** objects can be found in the following answers on StackOverflow: [[1]](http://stackoverflow.com/questions/38894044/print-beautiful-tables-for-h2o-models-in-r/39135080#39135080), [[2]](http://stackoverflow.com/questions/39397194/computing-p-values-in-spatial-econometric-models-why-are-there-inconsistencies/39479191#39479191), [[3]](http://stackoverflow.com/questions/36947477/how-can-i-use-texreg-1-36-4-for-a-relogit-model-estimated-using-zelig-v-5/36968738#36968738), [[4]](http://stackoverflow.com/questions/39143747/how-to-use-texreg-after-clmm-i-want-to-extract-random-effect-components/39507751#39507751), [[5]](http://stackoverflow.com/questions/40176607/r-how-to-get-a-proper-latex-regression-table-from-a-dataframe/40197961#40197961).

## Installation

The last stable release can be installed from CRAN:
``` r
install.packages("texreg")
```
To install the latest development version from GitHub, use the `remotes` package:
``` r
remotes::install_github("leifeld/texreg")
```

## Contribute to the project

Please feel free to report bugs or suggested enhancements using the [issue tracker](http://github.com/leifeld/texreg/issues) and propose solutions for known bugs using a pull request. Please observe and follow the code formatting in the texreg package when doing so, and add (or update) `testthat` unit tests to your pull requests if possible.

[![CRAN check and test](https://github.com/leifeld/texreg/actions/workflows/CRAN%20check%20and%20test.yaml/badge.svg)](https://github.com/leifeld/texreg/actions/workflows/CRAN%20check%20and%20test.yaml)
[![cran version](https://www.r-pkg.org/badges/version/texreg)](https://cran.r-project.org/package=texreg)
[![downloads](https://cranlogs.r-pkg.org/badges/texreg)](http://cranlogs.r-pkg.org/badges/texreg)
[![total downloads](https://cranlogs.r-pkg.org/badges/grand-total/texreg)](http://cranlogs.r-pkg.org/badges/grand-total/texreg)
[![Research software impact](http://depsy.org/api/package/cran/texreg/badge.svg)](http://depsy.org/package/r/texreg)
[![Coverage status](https://codecov.io/gh/leifeld/texreg/branch/master/graph/badge.svg)](https://codecov.io/github/leifeld/texreg?branch=master)
