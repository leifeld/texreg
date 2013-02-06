\name{texreg-package}
\docType{package}
\alias{texreg-package}
\title{Conversion of \R regression output to LaTeX tables.}
\description{texreg converts coefficients, standard errors, significance stars, 
    and goodness-of-fit statistics of statistical models into LaTeX tables. 
    Either a single model or a list of models can be combined in a single 
    table. The output is customizable. New model types can be easily 
    implemented.}
\details{
Several packages like apsrtable, xtable, outreg, stargazer, and memisc are 
available for typesetting R regression output as LaTeX tables. However, they 
cannot process certain model types (like ergm or mer objects), and some of 
these packages do not support multiple models within one single table. The 
texreg package, in contrast, accepts single models or lists of several models, 
and it can handle many different model types, including options for creating 
custom extensions for one's own model types. If several models are submitted, 
they are merged by row labels of the coefficients, and they are inserted into 
the final table as separate columns. Beside LaTeX output, texreg can also print 
nicely formatted regression tables to the R console screen or as HTML code. 
Nicely formatted tables can also be saved as MS Word documents.

To display citation information, execute \code{citation("texreg")}.
}
\seealso{
\code{\link{texreg}}
}
\author{
Philip Leifeld (\url{http://www.philipleifeld.de})
}
\keyword{Programming|programming|interface}
\keyword{Programming|IO}
