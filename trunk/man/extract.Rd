\name{extract}
\alias{extract}
\alias{extract.lm}
\alias{extract.glm}
\alias{extract.gls}
\alias{extract.lme}
\alias{extract.ergm}
\alias{extract.pmg}
\alias{extract.plm}
\alias{extract.rq}
\alias{extract.lrm}
\alias{extract.clogit}
\title{Extract coefficients and GOF measures from a statistical object}
\description{Extract coefficients and GOF measures from a statistical object.}
\usage{
extract(model, ...)
extract.lm(model)
extract.glm(model)
extract.gls(model)
extract.lme(model)
extract.ergm(model)
extract.pmg(model)
extract.plm(model)
extract.rq(model)
extract.lrm(model)
extract.clogit(model)
}
\details{
This is a generic function which extracts coefficients and GOF measures from 
statistical objects. There are several extract functions for the specific 
model types, which are called by the generic extract function if it encounters 
a model known to be handled by the specific function. The output is useful 
for the \link{texreg} function.
}
\arguments{
\item{model}{ A statistical model object. }
\item{...}{ Custom parameters. }
}
\seealso{
\code{\link{texreg-package} \link{texreg}}
}
\author{
Philip Leifeld (\url{http://www.philipleifeld.de})
}
