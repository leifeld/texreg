\name{extract}
\alias{extract}
\alias{extract.clm}
\alias{extract.clogit}
\alias{extract.coxph}
\alias{extract.coxph.penal}
\alias{extract.dynlm}
\alias{extract.ergm}
\alias{extract.gee}
\alias{extract.glm}
\alias{extract.glmerMod}
\alias{extract.gls}
\alias{extract.gmm}
\alias{extract.lm}
\alias{extract.lme}
\alias{extract.lmerMod}
\alias{extract.lmrob}
\alias{extract.lnam}
\alias{extract.lrm}
\alias{extract.mer}
\alias{extract.negbin}
\alias{extract.nlmerMod}
\alias{extract.plm}
\alias{extract.pmg}
\alias{extract.polr}
\alias{extract.Relogit}
\alias{extract.rlm}
\alias{extract.rq}
\alias{extract.sclm}
\alias{extract.simex}
\alias{extract.stergm}
\alias{extract.svyglm}
\alias{extract.systemfit}
\alias{extract.tobit}
\title{Extract coefficients and GOF measures from a statistical object}
\description{Extract coefficients and GOF measures from a statistical object.}
\usage{
extract(model, ...)
extract.clm(model, include.thresholds=TRUE, include.aic=TRUE, include.bic=TRUE, 
    include.loglik=TRUE, include.nobs=TRUE, ...)
extract.clogit(model, include.aic=TRUE, include.rsquared=TRUE, 
    include.maxrs=TRUE, include.events=TRUE, include.nobs=TRUE, 
    include.missings=TRUE, ...)
extract.coxph(model, include.aic=TRUE, include.rsquared=TRUE, 
    include.maxrs=TRUE, include.events=TRUE, include.nobs=TRUE, 
    include.missings=TRUE, include.zph=TRUE, ...)
extract.coxph.penal(model, include.aic=TRUE, include.rsquared=TRUE, 
    include.maxrs=TRUE, include.events=TRUE, include.nobs=TRUE, 
    include.missings=TRUE, include.zph=TRUE, ...)
extract.dynlm(model, include.rsquared=TRUE, include.adjrs=TRUE, 
    include.nobs=TRUE, ...)
extract.ergm(model, include.aic=TRUE, include.bic=TRUE, include.loglik=TRUE, 
    ...)
extract.gee(model, robust=TRUE, include.dispersion=TRUE, include.nobs=TRUE, ...)
extract.glm(model, include.aic=TRUE, include.bic=TRUE, include.loglik=TRUE, 
    include.deviance=TRUE, include.nobs=TRUE, ...)
extract.glmerMod(model, include.pvalues=FALSE, include.aic=TRUE, 
    include.bic=TRUE, include.loglik=TRUE, include.deviance=TRUE, 
    include.nobs=TRUE, include.groups=TRUE, include.variance=TRUE, 
    mcmc.pvalues=FALSE, mcmc.size=5000, ...)
extract.gls(model, include.aic=TRUE, include.bic=TRUE, include.loglik=TRUE, 
    include.nobs=TRUE, ...)
extract.gmm(model, include.obj.fcn=TRUE, include.overidentification=FALSE, 
    include.nobs=TRUE, ...)
extract.lm(model, include.rsquared=TRUE, include.adjrs=TRUE, include.nobs=TRUE, 
    ...)
extract.lme(model, include.aic=TRUE, include.bic=TRUE, include.loglik=TRUE, 
    include.nobs=TRUE, ...)
extract.lmerMod(model, include.pvalues=FALSE, include.aic=TRUE, 
    include.bic=TRUE, include.loglik=TRUE, include.deviance=TRUE, 
    include.nobs=TRUE, include.groups=TRUE, include.variance=TRUE, 
    mcmc.pvalues=FALSE, mcmc.size=5000, ...)
extract.lmrob(model, include.nobs=TRUE, ...)
extract.lnam(model, include.rsquared=TRUE, include.adjrs=TRUE, 
    include.aic=TRUE, include.bic=TRUE, include.loglik=TRUE, ...)
extract.lrm(model, include.pseudors=TRUE, include.lr=TRUE, include.nobs=TRUE, 
    ...)
extract.mer(model, include.pvalues=FALSE, include.aic=TRUE, 
    include.bic=TRUE, include.loglik=TRUE, include.deviance=TRUE, 
    include.nobs=TRUE, include.groups=TRUE, include.variance=TRUE, 
    mcmc.pvalues=FALSE, mcmc.size=5000, ...)
extract.negbin(model, include.aic=TRUE, include.bic=TRUE, include.loglik=TRUE, 
    include.deviance=TRUE, include.nobs=TRUE, ...)
extract.nlmerMod(model, include.pvalues=FALSE, include.aic=TRUE, 
    include.bic=TRUE, include.loglik=TRUE, include.deviance=TRUE, 
    include.nobs=TRUE, include.groups=TRUE, include.variance=TRUE, 
    mcmc.pvalues=FALSE, mcmc.size=5000, ...)
extract.plm(model, include.rsquared=TRUE, include.adjrs=TRUE, 
    include.nobs=TRUE, ...)
extract.pmg(model, include.nobs=TRUE, ...)
extract.polr(model, include.thresholds=FALSE, include.aic=TRUE, 
    include.bic=TRUE, include.loglik=TRUE, include.deviance=TRUE, 
    include.nobs=TRUE, ...)
extract.Relogit(model, include.aic=TRUE, include.bic=TRUE, include.loglik=TRUE, 
    include.deviance=TRUE, include.nobs=TRUE, ...)
extract.rlm(model, include.nobs=TRUE, ...)
extract.rq(model, include.nobs=TRUE, include.percentile=TRUE, ...)
extract.sclm(model, include.thresholds=TRUE, include.aic=TRUE, include.bic=TRUE,
    include.loglik=TRUE, include.nobs=TRUE, ...)
extract.simex(model, jackknife=TRUE, include.nobs=TRUE, ...)
extract.stergm(model, beside=FALSE, include.formation=TRUE, 
    include.dissolution=TRUE, include.nvertices=TRUE, include.aic=FALSE, 
    include.bic=FALSE, include.loglik=FALSE, ...)
extract.svyglm(model, include.aic=FALSE, include.bic=FALSE, 
    include.loglik=FALSE,include.deviance=TRUE, include.dispersion=TRUE, 
    include.nobs=TRUE, ...)
extract.systemfit(model, include.rsquared=TRUE, include.adjrs=TRUE, 
    include.nobs=TRUE, ...)
extract.tobit(model, include.aic=TRUE, include.bic=TRUE, include.loglik=TRUE, 
    include.deviance=TRUE, include.nobs=FALSE, include.censnobs=TRUE, 
    include.wald=TRUE, ...)
}
\details{
extract() is a generic function which extracts coefficients and GOF measures from 
statistical objects. There are several extract functions for the specific 
model types, which are called by the generic extract function if it encounters 
a model known to be handled by the specific function. The output is used by 
the \link{texreg} function.

The various extract functions can also be used directly on a statistical model 
to convert them into texreg objects.
}
\arguments{
\item{model}{ A statistical model object. }
\item{beside}{ If available: should the model terms be arranged below each other or beside each other (default)? For example, in a stergm model, the formation and dissolution coefficients can be arranged in two different ways. }
\item{include.adjrs}{If available: should the adjusted R-squared be reported? }
\item{include.aic}{ If available: should Akaike's information criterion (AIC) be reported? }
\item{include.bic}{ If available: should the Bayesian information criterion (BIC) be reported? }
\item{include.censnobs}{ If available: should the total, right-censored, left-censored, and uncensored number of observations be reported? }
\item{include.deviance}{ If available: should the deviance be reported? }
\item{include.dispersion}{ If available: should the dispersion or scale parameter be reported? }
\item{include.dissolution}{ If available: should the coefficients for the dissolution phase in a STERGM be reported? }
\item{include.events}{ If available: should the number of events be reported (in survival models)? }
\item{include.formation}{ If available: should the coefficients for the formation phase in a STERGM be reported? }
\item{include.groups}{ If available: should the number of groups be reported? }
\item{include.loglik}{ If available: should the log-likelihood be reported? }
\item{include.lr}{ If available: should the likelihood ratio test be reported? }
\item{include.maxrs}{If available: should the maximum possible R-squared be reported? }
\item{include.missings}{ If available: should the number of missing observations be reported (in survival models)? }
\item{include.nobs}{ If available: should the number of observations be reported? }
\item{include.nvertices}{ If available: should the number of vertices be reported in a statistical network model? }
\item{include.obj.fcn}{ If available: should the value of the objective function (= criterion function) be reported (for \code{gmm} objects)? More precisely, this returns \code{E(g)var(g)^{-1}E(g)}. }
\item{include.overidentification}{ If available: should the J-test for overidentification be reported (for \code{gmm} objects)? }
\item{include.percentile}{ If available: should the percentile (tau) be reported? }
\item{include.pseudors}{ If available: should the pseudo R-squared be reported? }
\item{include.pvalues}{ If available: should the p-values be reported (not recommended for lme4 models)? }
\item{include.rsquared}{ If available: should R-squared be reported? }
\item{include.thresholds}{ If available: should the threshold parameters (that is, the intercepts for the class boundaries) be reported in ordinal models? }
\item{include.variance}{ If available: should group variances be reported? }
\item{include.wald}{ If available: should the Wald statistic be included? }
\item{include.zph}{ If available: should the Cox proportional hazards assumption be tested (resulting in a p value indicating whether the proportional hazards assumption of the model is violated)? }
\item{jackknife}{ If available: use Jackknife variance instead of Asymptotic variance. }
\item{mcmc.pvalues}{ In linear mixed effects models: compute MCMC-based p-values instead of naive p-values. This will only affect the output if the argument \code{include.pvalues=TRUE} is also set. Warning: computing MCMC-based p-values may take some time. }
\item{mcmc.size}{ In linear mixed effects models: the MCMC sample size from which p-values are derived (if the arguments \code{include.pvalues=TRUE} and \code{mcmc.pvalues=TRUE} are also set). Note: high values may take considerable computing time. }
\item{robust}{ If available: report robust instead of native standard errors. }
\item{...}{ Custom parameters. }
}
\seealso{
\code{\link{texreg-package} \link{texreg} \link{extract-methods}}
}
\author{
Philip Leifeld (\url{http://www.philipleifeld.de})
}
