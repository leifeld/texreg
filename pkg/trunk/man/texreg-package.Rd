\name{texreg-package}
\docType{package}
\alias{texreg-package}
\title{Conversion of \R regression output to LaTeX tables.}
\description{texreg converts coefficients, standard errors, significance stars, 
    and goodness-of-fit statistics of lme objects (that is, linear mixed 
    effects models) and ergm objects (that is, exponential random graph models) 
    into LaTeX tables. Either a single model or a list of models can be 
    combined in a single table. The output is customizable.}
\details{
Several packages like xtable, outreg and memsic are available for typesetting 
R regression output as LaTeX tables. However, these packages do not work with 
ergm objects (that is, exponential random graph models) and lme objects (that 
is, linear mixed effects models). The texreg package, in contrast, is able to 
convert lme and ergm output to LaTeX tables. It accepts single models or lists 
of several models. If several models are submitted, they are merged by row 
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
