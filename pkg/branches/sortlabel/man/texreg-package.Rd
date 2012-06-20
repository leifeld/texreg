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
Several packages like xtable, outreg and memsic are available for typesetting 
R regression output as LaTeX tables. However, they cannot process certain model 
types (like ergm or lme objects), and some of these packages do not support 
multiple models within one single table. The texreg package, in contrast, 
accepts single models or lists of several models, and it can handle ergm and 
lme objects. If several models are submitted, they are merged by row 
labels of the coefficients, and they are inserted into the final table as 
separate columns.

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
