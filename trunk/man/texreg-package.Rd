\name{texreg-package}
\docType{package}
\alias{texreg-package}
\title{Conversion of statistical model output in \R to LaTeX and HTML tables}
\description{
texreg converts coefficients, standard errors, significance stars, and 
goodness-of-fit statistics of statistical models into LaTeX tables or HTML 
tables/MS Word documents or to nicely formatted screen output for the R console 
for easy model comparison. A list of several models can be combined in a single 
table. The output is highly customizable. New model types can be easily 
implemented.}
\details{
Several packages like apsrtable, memisc, outreg, stargazer, or xtable are 
available for typesetting \R regression output as LaTeX tables. However, 
texreg supports more flexible handling of new model types, supports multiple 
models side-by-side, supports confidence intervals and standard errors alike,
has many options for customization, and beside LaTeX output, it can also 
print tables to the \R console screen, save them as HTML or MS Word documents, 
or plot them as coefficient plots. If several models are submitted, they are 
merged by row labels of the coefficients, and they are inserted into 
the final table as separate columns. The package works with report generation 
tools like Sweave or knitr (including options for LaTeX, HTML, and Markdown). 
To display citation information, execute \code{citation("texreg")}.
}
\seealso{
\code{\link{texreg} \link{plotreg}}
}
\author{
Philip Leifeld (\url{http://www.philipleifeld.de})
}
\keyword{print}
\keyword{misc}
\keyword{utilities}
\keyword{IO}
\keyword{programming|interface}
