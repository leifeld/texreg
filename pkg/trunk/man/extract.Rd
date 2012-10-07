\name{extract}
\alias{extract}
\alias{extract.clogit}
\alias{extract.ergm}
\alias{extract.glm}
\alias{extract.gls}
\alias{extract.lm}
\alias{extract.lme}
\alias{extract.lnam}
\alias{extract.lrm}
\alias{extract.plm}
\alias{extract.pmg}
\alias{extract.rq}
\alias{extract.systemfit}
\title{Extract coefficients and GOF measures from a statistical object}
\description{Extract coefficients and GOF measures from a statistical object.}
\usage{
extract(model, ...)
extract.clogit(model)
extract.ergm(model)
extract.glm(model)
extract.gls(model)
extract.lm(model)
extract.lme(model)
extract.lnam{model}
extract.lrm(model)
extract.plm(model)
extract.pmg(model)
extract.rq(model)
extract.systemfit(model)
}
\details{
This is a generic function which extracts coefficients and GOF measures from 
statistical objects. There are several extract functions for the specific 
model types, which are called by the generic extract function if it encounters 
a model known to be handled by the specific function. The output is used by 
the \link{texreg} function.
}
\arguments{
\item{model}{ A statistical model object. }
\item{...}{ Custom parameters. }
}
\seealso{
\code{\link{texreg-package} \link{texreg} \link{extract-methods}}
}
\author{
Philip Leifeld (\url{http://www.philipleifeld.de})
}
