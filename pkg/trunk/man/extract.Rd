\name{extract}
\alias{extract}
\alias{extract.aftreg}
\alias{extract.betareg}
\alias{extract.clm}
\alias{extract.clogit}
\alias{extract.coxph}
\alias{extract.coxph.penal}
\alias{extract.dynlm}
\alias{extract.ergm}
\alias{extract.gam}
\alias{extract.gee}
\alias{extract.glm}
\alias{extract.glmerMod}
\alias{extract.gls}
\alias{extract.gmm}
\alias{extract.hurdle}
\alias{extract.ivreg}
\alias{extract.lm}
\alias{extract.lme}
\alias{extract.lmerMod}
\alias{extract.lmrob}
\alias{extract.lnam}
\alias{extract.lrm}
\alias{extract.maBina}
\alias{extract.mer}
\alias{extract.multinom}
\alias{extract.negbin}
\alias{extract.nlmerMod}
\alias{extract.phreg}
\alias{extract.plm}
\alias{extract.pmg}
\alias{extract.polr}
\alias{extract.Relogit}
\alias{extract.rem.dyad}
\alias{extract.rlm}
\alias{extract.rq}
\alias{extract.sclm}
\alias{extract.simex}
\alias{extract.stergm}
\alias{extract.survreg}
\alias{extract.survreg.penal}
\alias{extract.svyglm}
\alias{extract.systemfit}
\alias{extract.tobit}
\alias{extract.weibreg}
\alias{extract.zeroinfl}
\title{Extract coefficients and GOF measures from a statistical object}
\description{Extract coefficients and GOF measures from a statistical object.}
\usage{
extract(model, ...)

extract.aftreg(model, include.loglik = TRUE, include.lr = TRUE, 
    include.nobs = TRUE, include.events = TRUE, include.trisk = TRUE, 
    ...)

extract.betareg(model, include.precision = TRUE, 
    include.pseudors = TRUE, include.loglik = TRUE, 
    include.nobs = TRUE, ...)

extract.clm(model, include.thresholds = TRUE, include.aic = TRUE, 
    include.bic=TRUE, include.loglik = TRUE, include.nobs = TRUE, 
    ...)

extract.clogit(model, include.aic = TRUE, include.rsquared = TRUE, 
    include.maxrs = TRUE, include.events = TRUE, 
    include.nobs = TRUE, include.missings = TRUE, ...)

extract.coxph(model, include.aic = TRUE, include.rsquared = TRUE, 
    include.maxrs=TRUE, include.events = TRUE, 
    include.nobs = TRUE, include.missings = TRUE, 
    include.zph = TRUE, ...)

extract.coxph.penal(model, include.aic = TRUE, 
    include.rsquared = TRUE, include.maxrs = TRUE, 
    include.events = TRUE, include.nobs = TRUE, 
    include.missings = TRUE, include.zph = TRUE, ...)

extract.dynlm(model, include.rsquared = TRUE, include.adjrs = TRUE, 
    include.nobs = TRUE, ...)

extract.ergm(model, include.aic = TRUE, include.bic = TRUE, 
    include.loglik=TRUE, ...)

extract.gam(model, include.smooth = TRUE, include.aic = TRUE, 
    include.bic = TRUE, include.loglik = TRUE, 
    include.deviance = TRUE, include.dev.expl = TRUE, 
    include.dispersion = TRUE, include.rsquared = TRUE, 
    include.gcv = TRUE, include.nobs = TRUE, 
    include.nsmooth = TRUE, ...)

extract.gee(model, robust = TRUE, include.dispersion = TRUE, 
    include.nobs = TRUE, ...)

extract.glm(model, include.aic = TRUE, include.bic = TRUE, 
    include.loglik = TRUE, include.deviance = TRUE, 
    include.nobs = TRUE, ...)

extract.glmerMod(model, include.pvalues = FALSE, include.aic = TRUE, 
    include.bic = TRUE, include.loglik = TRUE, 
    include.deviance = TRUE, include.nobs = TRUE, 
    include.groups = TRUE, include.variance = TRUE, 
    mcmc.pvalues = FALSE, mcmc.size = 5000, ...)

extract.gls(model, include.aic = TRUE, include.bic = TRUE, 
    include.loglik = TRUE, include.nobs = TRUE, ...)

extract.gmm(model, include.obj.fcn = TRUE, 
    include.overidentification = FALSE, include.nobs = TRUE, ...)

extract.hurdle(model, beside = FALSE, include.count = TRUE, 
    include.zero = TRUE, include.aic = TRUE, include.loglik = TRUE, 
    include.nobs = TRUE, ...)

extract.ivreg(model, include.rsquared = TRUE, include.adjrs = TRUE, 
    include.nobs = TRUE, ...)

extract.lm(model, include.rsquared = TRUE, include.adjrs = TRUE, 
    include.nobs = TRUE, ...)

extract.lme(model, include.aic = TRUE, include.bic = TRUE, 
    include.loglik = TRUE, include.nobs = TRUE, ...)

extract.lmerMod(model, include.pvalues = FALSE, include.aic = TRUE, 
    include.bic = TRUE, include.loglik = TRUE, 
    include.deviance = TRUE, include.nobs = TRUE, 
    include.groups = TRUE, include.variance = TRUE, 
    mcmc.pvalues = FALSE, mcmc.size = 5000, ...)

extract.lmrob(model, include.nobs = TRUE, ...)

extract.lnam(model, include.rsquared = TRUE, include.adjrs = TRUE, 
    include.aic = TRUE, include.bic = TRUE, include.loglik = TRUE, 
    ...)

extract.lrm(model, include.pseudors = TRUE, include.lr = TRUE, 
    include.nobs = TRUE, ...)

extract.maBina(model, ...)

extract.mer(model, include.pvalues = FALSE, include.aic = TRUE, 
    include.bic = TRUE, include.loglik = TRUE, 
    include.deviance = TRUE, include.nobs = TRUE, 
    include.groups = TRUE, include.variance = TRUE, 
    mcmc.pvalues = FALSE, mcmc.size = 5000, ...)

extract.multinom(model, include.pvalues = TRUE, include.aic = TRUE, 
    include.bic = TRUE, include.loglik = TRUE, 
    include.deviance = TRUE, include.nobs = TRUE, ...)

extract.negbin(model, include.aic = TRUE, include.bic = TRUE, 
    include.loglik = TRUE, include.deviance = TRUE, 
    include.nobs = TRUE, ...)

extract.nlmerMod(model, include.pvalues = FALSE, include.aic = TRUE, 
    include.bic = TRUE, include.loglik = TRUE, 
    include.deviance = TRUE, include.nobs = TRUE, 
    include.groups = TRUE, include.variance = TRUE, 
    mcmc.pvalues=FALSE, mcmc.size=5000, ...)

extract.phreg(model, include.loglik = TRUE, include.lr = TRUE, 
    include.nobs = TRUE, include.events = TRUE, include.trisk = TRUE, 
    ...)

extract.plm(model, include.rsquared = TRUE, include.adjrs = TRUE, 
    include.nobs = TRUE, ...)

extract.pmg(model, include.nobs = TRUE, ...)

extract.polr(model, include.thresholds = FALSE, include.aic = TRUE, 
    include.bic = TRUE, include.loglik = TRUE, 
    include.deviance = TRUE, include.nobs = TRUE, ...)

extract.Relogit(model, include.aic = TRUE, include.bic = TRUE, 
    include.loglik = TRUE, include.deviance = TRUE, 
    include.nobs = TRUE, ...)

extract.rem.dyad(model, include.nvertices = TRUE, 
    include.events = TRUE, include.aic = TRUE, include.aicc = TRUE, 
    include.bic = TRUE, ...)

extract.rlm(model, include.nobs = TRUE, ...)

extract.rq(model, include.nobs = TRUE, include.percentile = TRUE, 
    ...)

extract.sclm(model, include.thresholds = TRUE, include.aic = TRUE, 
    include.bic = TRUE, include.loglik = TRUE, include.nobs = TRUE, 
    ...)

extract.simex(model, jackknife = TRUE, include.nobs = TRUE, ...)

extract.stergm(model, beside = FALSE, include.formation = TRUE, 
    include.dissolution = TRUE, include.nvertices = TRUE, 
    include.aic = FALSE, include.bic = FALSE, 
    include.loglik = FALSE, ...)

extract.survreg(model, include.aic = TRUE, include.bic = TRUE, 
    include.loglik = TRUE, include.deviance = TRUE, 
    include.nobs = TRUE, ...)

extract.survreg.penal(model, include.aic = TRUE, 
    include.bic = TRUE, include.loglik = TRUE, 
    include.deviance = TRUE, include.nobs = TRUE, ...)

extract.svyglm(model, include.aic = FALSE, include.bic = FALSE, 
    include.loglik = FALSE, include.deviance = TRUE, 
    include.dispersion = TRUE, include.nobs = TRUE, ...)

extract.systemfit(model, include.rsquared = TRUE, 
    include.adjrs = TRUE, include.nobs = TRUE, ...)

extract.tobit(model, include.aic = TRUE, include.bic = TRUE, 
    include.loglik = TRUE, include.deviance = TRUE, 
    include.nobs = FALSE, include.censnobs = TRUE, include.wald=TRUE, 
    ...)

extract.weibreg(model, include.loglik = TRUE, include.lr = TRUE, 
    include.nobs = TRUE, include.events = TRUE, include.trisk = TRUE, 
    ...)

extract.zeroinfl(model, beside = FALSE, include.count = TRUE, 
    include.zero = TRUE, include.aic = TRUE, include.loglik = TRUE, 
    include.nobs = TRUE, ...)
}
\details{
extract() is a generic function which extracts coefficients and GOF measures 
from statistical objects. There are several extract functions for the specific 
model types, which are called by the generic extract function if it encounters 
a model known to be handled by the specific function. The output is used by 
the \link{texreg} function.

The various extract functions can also be used directly on a statistical model 
to convert them into texreg objects.
}
\arguments{
\item{model}{ A statistical model object. }
\item{beside}{ If available: should the model terms be arranged below each other or beside each other (default)? For example, in a stergm model, the formation and dissolution coefficients can be arranged in two columns of the table. }
\item{include.adjrs}{If available: should the adjusted R-squared be reported? }
\item{include.aic}{ If available: should Akaike's information criterion (AIC) be reported? }
\item{include.aicc}{ If available: should AICc be reported? This is a version of AIC with a correction for finite sample sizes. }
\item{include.bic}{ If available: should the Bayesian information criterion (BIC) be reported? }
\item{include.censnobs}{ If available: should the total, right-censored, left-censored, and uncensored number of observations be reported? }
\item{include.count}{ If available: should the count model of a zero-inflated or hurdle regression be included in the coefficients block (before the zero-inflation or zero hurdle model)? }
\item{include.dev.expl}{ If available: should the deviance explained be reported? }
\item{include.deviance}{ If available: should the deviance be reported? }
\item{include.dispersion}{ If available: should the dispersion or scale parameter be reported? }
\item{include.dissolution}{ If available: should the coefficients for the dissolution phase in a STERGM be reported? }
\item{include.events}{ If available: should the number of events be reported (in survival models)? }
\item{include.formation}{ If available: should the coefficients for the formation phase in a STERGM be reported? }
\item{include.gcv}{ If available: should the GCV score be reported (in GAMs)? }
\item{include.groups}{ If available: should the number of groups be reported? }
\item{include.loglik}{ If available: should the log-likelihood be reported? }
\item{include.lr}{ If available: should the likelihood ratio test be reported? }
\item{include.maxrs}{If available: should the maximum possible R-squared be reported? }
\item{include.missings}{ If available: should the number of missing observations be reported (in survival models)? }
\item{include.nobs}{ If available: should the number of observations be reported? }
\item{include.nsmooth}{ If available: should the number of smooth terms be reported (in GAMs)? }
\item{include.nvertices}{ If available: should the number of vertices be reported in a statistical network model? }
\item{include.obj.fcn}{ If available: should the value of the objective function (= criterion function) be reported (for \code{gmm} objects)? More precisely, this returns \code{E(g)var(g)^{-1}E(g)}. }
\item{include.overidentification}{ If available: should the J-test for overidentification be reported (for \code{gmm} objects)? }
\item{include.percentile}{ If available: should the percentile (tau) be reported? }
\item{include.precision}{ If available: should the precision estimates of a betareg fit (the phi coefficients) be reported as part of the coefficients block? }
\item{include.pseudors}{ If available: should the pseudo R-squared be reported? }
\item{include.pvalues}{ If available: should the p values be reported (naive p values are not recommended for lme4 models, but see also the \code{mcmc.pvalues} argument)? }
\item{include.rsquared}{ If available: should R-squared be reported? }
\item{include.smooth}{ If available: should the smooth terms of a GAM be reported? If they are reported, the EDF value is reported as the coefficient, and DF is included in parentheses (not standard errors because a chi-square test is used for the smooth terms). }
\item{include.thresholds}{ If available: should the threshold parameters (that is, the intercepts for the class boundaries) be reported in ordinal models? }
\item{include.trisk}{ If available: should the total time at risk be reported (in event-history models)? }
\item{include.variance}{ If available: should group variances be reported? }
\item{include.wald}{ If available: should the Wald statistic be included? }
\item{include.zero}{ If available: should the zero-inflation model of a zero-inflated regression or the zero hurdle model of a hurdle regression be included in the coefficients block (after the count model)? }
\item{include.zph}{ If available: should the Cox proportional hazards assumption be tested (resulting in a p value indicating whether the proportional hazards assumption of the model is violated)? }
\item{jackknife}{ If available: use Jackknife variance instead of Asymptotic variance. }
\item{mcmc.pvalues}{ In linear mixed effects models: compute MCMC-based p values instead of naive p values. This will only affect the output if the argument \code{include.pvalues=TRUE} is also set. Warning: computing MCMC-based p values may take some time. }
\item{mcmc.size}{ In linear mixed effects models: the MCMC sample size from which p values are derived (if the arguments \code{include.pvalues = TRUE} and \code{mcmc.pvalues = TRUE} are also set). Note: high values may take considerable computing time. }
\item{robust}{ If available: report robust instead of native standard errors. }
\item{...}{ Custom parameters. }
}
\seealso{
\code{\link{texreg-package} \link{texreg} \link{extract-methods}}
}
\author{
Philip Leifeld (\url{http://www.philipleifeld.de})
}
