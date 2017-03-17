# texreg

Conversion of R Regression Output to LaTeX or HTML Tables.

Converts coefficients, standard errors, significance stars, and goodness-of-fit statistics of statistical models into LaTeX tables or HTML tables/MS Word documents or to nicely formatted screen output for the R console for easy model comparison. A list of several models can be combined in a single table. The output is highly customizable. New model types can be easily implemented.

[![Build Status](https://travis-ci.org/leifeld/texreg.svg?branch=master)](https://travis-ci.org/leifeld/texreg)
[![cran version](http://www.r-pkg.org/badges/version/texreg)](https://cran.r-project.org/package=texreg)
[![downloads](http://cranlogs.r-pkg.org/badges/texreg)](http://cranlogs.r-pkg.org/badges/texreg)
[![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/texreg)](http://cranlogs.r-pkg.org/badges/grand-total/texreg)
[![Research software impact](http://depsy.org/api/package/cran/texreg/badge.svg)](http://depsy.org/package/r/texreg)

## Documentation

Details on **texreg** can be found in the following article:

Leifeld, Philip (2013): **texreg**: Conversion of Statistical Model Output in R to LaTeX and HTML Tables. Journal of Statistical Software 55(8): 1-24. doi:[10.18637/jss.v055.i08](http://dx.doi.org/10.18637/jss.v055.i08)

An updated version of this paper is included as a [package vignette](https://cran.r-project.org/web/packages/texreg/vignettes/texreg.pdf). Additional details on how to write custom extensions for **texreg** and manipulate **texreg** objects can be found in the following answers on StackOverflow: [[1]](http://stackoverflow.com/questions/38894044/print-beautiful-tables-for-h2o-models-in-r/39135080#39135080), [[2]](http://stackoverflow.com/questions/39397194/computing-p-values-in-spatial-econometric-models-why-are-there-inconsistencies/39479191#39479191), [[3]](http://stackoverflow.com/questions/36947477/how-can-i-use-texreg-1-36-4-for-a-relogit-model-estimated-using-zelig-v-5/36968738#36968738), [[4]](http://stackoverflow.com/questions/39143747/how-to-use-texreg-after-clmm-i-want-to-extract-random-effect-components/39507751#39507751), [[5]](http://stackoverflow.com/questions/40176607/r-how-to-get-a-proper-latex-regression-table-from-a-dataframe/40197961#40197961).


## Bugs or feature requests

Please feel free to report bugs or possible enhancements using the [issue tracker](http://github.com/leifeld/texreg/issues) and propose solutions for known bugs using a pull request. I am overwhelmed with other work at the moment and need your help in fixing bugs or adding new functions! Important: I will only consider pull requests that address an issue described in the issue tracker, so please create an issue first before fixing it.
