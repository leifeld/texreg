\name{createTexreg}
\alias{createTexreg}
\title{Create a texreg object}
\description{Create a texreg object with coefficients and GOF statistics.}
\usage{
createTexreg(coef.names, coef, se, pvalues = numeric(0), ci.low = numeric(0), 
    ci.up = numeric(0), gof.names = character(0), gof = numeric(0), 
    gof.decimal = logical(0))
}
\details{
This function creates a texreg object. A texreg object contains 
information about coefficients, standard errors, p values 
(optional), and about goodness-of-fit statistics. Instead of standard
errors and p values, a texreg object may also contain upper and
lower bounds of a confidence interval. texreg objects 
are used by the texreg command to create LaTeX tables and other
representations of the model results.
}
\arguments{
\item{coef.names}{ A vector of coefficient names. }
\item{coef}{ The coefficient values. }
\item{se}{ The standard errors. This is optional if the \code{ci.low} and \code{ci.up} slots are filled. }
\item{pvalues}{ The p-values of the model. This is optional. }
\item{ci.low}{ Lower bound of confidence interval (the actual values, not the confidence level). This is optional as long as \code{se} is available, but if it is provided, \code{ci.up} must also be provided. }
\item{ci.up}{ Upper bound of confidence interval (the actual values, not the confidence level). This is optional as long as \code{se} is available, but if it is provided, \code{ci.low} must also be provided. }
\item{gof.names}{ A vector of names of the goodness-of-fit statistics. }
\item{gof}{ A vector of goodness-of-fit statistics. }
\item{gof.decimal}{ A vector of boolean/logical values indicating for each GOF statistic if decimal places shall be used. This is optional. }
}
\seealso{
\code{\link{texreg-package} \link{texreg}}
}
\author{
Philip Leifeld (\url{http://www.philipleifeld.de})
}
\examples{
library(nlme)  #load library for fitting linear mixed effects models
model <- lme(distance ~ age, data = Orthodont, random = ~ 1)  #estimate model
coefficient.names <- rownames(summary(model)$tTable)  #extract coefficient names
coefficients <- summary(model)$tTable[, 1]  #extract coefficient values
standard.errors <- summary(model)$tTable[, 2]  #extract standard errors
significance <- summary(model)$tTable[, 5]  #extract p values

lik <- summary(model)$logLik  #extract log likelihood
aic <- summary(model)$AIC  #extract AIC
bic <- summary(model)$BIC  #extract BIC
n <- nobs(model)  #extract number of observations
gof <- c(aic, bic, lik, n)  #create a vector of GOF statistics
gof.names <- c("AIC", "BIC", "Log Likelihood", "Num. obs.")  #names of GOFs
decimal.places <- c(TRUE, TRUE, TRUE, FALSE)  #the last one is a count variable

#create the texreg object
tr <- createTexreg(
  coef.names = coefficient.names, 
  coef = coefficients, 
  se = standard.errors, 
  pvalues = significance, 
  gof.names = gof.names, 
  gof = gof, 
  gof.decimal = decimal.places
)
}
