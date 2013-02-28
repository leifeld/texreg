# The texreg package was written by Philip Leifeld.
# Please use the forum at http://r-forge.r-project.org/projects/texreg/ 
# for bug reports, help or feature requests.


# generic extract function
setGeneric("extract", function(model, ...) standardGeneric("extract"), 
    package="texreg")


# extension for clm objects
extract.clm <- function(model, include.thresholds=TRUE, include.aic=TRUE, 
    include.bic=TRUE, include.loglik=TRUE, include.nobs=TRUE, ...) {
  s <- summary(model, ...)
  
  tab <- s$coefficients
  thresh <- tab[rownames(tab) %in% names(s$aliased$alpha),]
  threshold.names <- rownames(thresh)
  threshold.coef <- thresh[,1]
  threshold.se <- thresh[,2]
  threshold.pval <- thresh[,4]
  beta <- tab[rownames(tab) %in% names(s$aliased$beta),]
  beta.names <- rownames(beta)
  beta.coef <- beta[,1]
  beta.se <- beta[,2]
  beta.pval <- beta[,4]
  if (include.thresholds==TRUE) {
    names <- c(beta.names, threshold.names)
    coef <- c(beta.coef, threshold.coef)
    se <- c(beta.se, threshold.se)
    pval <- c(beta.pval, threshold.pval)
  } else {
    names <- beta.names
    coef <- beta.coef
    se <- beta.se
    pval <- beta.pval
  }
  
  n <- nobs(model)
  lik <- logLik(model)[1]
  aic <- AIC(model)
  bic <- BIC(model)
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic==TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic==TRUE) {
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik==TRUE) {
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names=names, 
      coef=coef, 
      se=se, 
      pvalues=pval, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("clm", "ordinal"), 
    definition = extract.clm)

extract.sclm <- extract.clm
setMethod("extract", signature=className("sclm", "ordinal"), 
    definition = extract.clm)


# extension for coxph objects (survival package)
extract.coxph <- function(model, include.aic=TRUE, include.rsquared=TRUE, 
    include.maxrs=TRUE, include.events=TRUE, include.nobs=TRUE, 
    include.missings=TRUE, include.zph=TRUE, ...) {
  s <- summary(model, ...)
  
  coefficient.names <- rownames(s$coef)
  coefficients <- s$coef[,1]
  standard.errors <- s$coef[,3]
  significance <- s$coef[,5]
  
  aic <- extractAIC(model)[2]
  event <- model$nevent
  n <- model$n
  mis <- length(model$na.action)
  rs <- s$rsq[1]
  maxrs <- s$rsq[2]
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic==TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.rsquared==TRUE) {
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.maxrs==TRUE) {
    gof <- c(gof, maxrs)
    gof.names <- c(gof.names, "Max.\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.events==TRUE) {
    gof <- c(gof, event)
    gof.names <- c(gof.names, "Num.\ events")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.missings==TRUE) {
    gof <- c(gof, mis)
    gof.names <- c(gof.names, "Missings")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.zph==TRUE) {
    zph <- cox.zph(model)$table
    zph <- zph[length(zph[,1]), length(zph[1,])]
    gof <- c(gof, zph)
    gof.names <- c(gof.names, "PH test")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  
  tr <- createTexreg(
      coef.names=coefficient.names, 
      coef=coefficients, 
      se=standard.errors, 
      pvalues=significance, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("coxph", "survival"), 
    definition = extract.coxph)


# extension for coxph.penal objects (survival package)
extract.coxph.penal <- function(model, include.aic=TRUE, include.rsquared=TRUE,
    include.maxrs=TRUE, include.events=TRUE, include.nobs=TRUE, 
    include.missings=TRUE, include.zph=TRUE, ...) {
  
  coefficients <- coef(model, ...)
  coefficient.names <- names(coefficients)
  if (!is.null(model$naive.var)) {
    standard.errors <- sqrt(diag(model$naive.var))
  } else {
    standard.errors <- sqrt(diag(model$var))
  }
  significance <- 1 - pchisq((coefficients/standard.errors)^2, 1)

  aic <- extractAIC(model)[2]
  event <- model$nevent
  n <- model$n
  mis <- length(model$na.action)
  logtest <- -2 * (model$loglik[1] - model$loglik[2])
  rs <- 1 - exp( - logtest / model$n)
  maxrs <- 1 - exp((2 * model$loglik[1]) / model$n)

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic==TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.rsquared==TRUE) {
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.maxrs==TRUE) {
    gof <- c(gof, maxrs)
    gof.names <- c(gof.names, "Max.\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.events==TRUE) {
    gof <- c(gof, event)
    gof.names <- c(gof.names, "Num.\ events")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.missings==TRUE) {
    gof <- c(gof, mis)
    gof.names <- c(gof.names, "Missings")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.zph==TRUE) {
    zph <- cox.zph(model)$table
    zph <- zph[length(zph[,1]), length(zph[1,])]
    gof <- c(gof, zph)
    gof.names <- c(gof.names, "PH test")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  
  tr <- createTexreg(
      coef.names=coefficient.names, 
      coef=coefficients, 
      se=standard.errors, 
      pvalues=significance, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("coxph.penal", "survival"), 
    definition = extract.coxph.penal)


# extension for clogit objects (survival package)
extract.clogit <- function(model, include.aic=TRUE, include.rsquared=TRUE, 
    include.maxrs=TRUE, include.events=TRUE, include.nobs=TRUE, 
    include.missings=TRUE, ...) {
  s <- summary(model, ...)
  
  coefficient.names <- rownames(s$coef)
  coefficients <- s$coef[,1]
  standard.errors <- s$coef[,3]
  significance <- s$coef[,5]
  
  aic <- extractAIC(model)[2]
  event <- model$nevent
  n <- model$n
  mis <- length(model$na.action)
  rs <- s$rsq[1]
  maxrs <- s$rsq[2]
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic==TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.rsquared==TRUE) {
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.maxrs==TRUE) {
    gof <- c(gof, maxrs)
    gof.names <- c(gof.names, "Max.\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.events==TRUE) {
    gof <- c(gof, event)
    gof.names <- c(gof.names, "Num.\ events")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.missings==TRUE) {
    gof <- c(gof, mis)
    gof.names <- c(gof.names, "Missings")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names=coefficient.names, 
      coef=coefficients, 
      se=standard.errors, 
      pvalues=significance, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("clogit", "survival"), 
    definition = extract.clogit)


# extension for ergm objects
extract.ergm <- function(model, include.aic=TRUE, include.bic=TRUE, 
    include.loglik=TRUE, ...) {
  s <- summary(model, ...)
  
  coefficient.names <- rownames(s$coefs)
  coefficients <- s$coefs[,1]
  standard.errors <- s$coefs[,2]
  significance <- s$coefs[,4]
  
  lik <- model$mle.lik[1]
  aic <- s$aic
  bic <- s$bic
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic==TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic==TRUE) {
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik==TRUE) {
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  
  tr <- createTexreg(
      coef.names=coefficient.names, 
      coef=coefficients, 
      se=standard.errors, 
      pvalues=significance, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("ergm", "ergm"), 
    definition = extract.ergm)


# extension for gee objects (gee package)
extract.gee <- function(model, robust=TRUE, include.dispersion=TRUE, 
    include.nobs=TRUE, ...) {
  s <- summary(model, ...)
  
  names <- rownames(coef(s))
  co <- coef(s)[,1]
  if (robust==TRUE) {
    se <- coef(s)[,4]
    zval <- coef(s)[,5]
  } else {
    se <- coef(s)[,2]
    zval <- coef(s)[,3]
  }
  pval <- pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  
  n <- nobs(model)
  disp <- s$scale
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.dispersion==TRUE) {
    gof <- c(gof, disp)
    gof.names <- c(gof.names, "Dispersion")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names=names, 
      coef=co, 
      se=se, 
      pvalues=pval, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("gee", "gee"), 
    definition = extract.gee)


# extension for glm objects
extract.glm <- function(model, include.aic=TRUE, include.bic=TRUE, 
    include.loglik=TRUE, include.deviance=TRUE, include.nobs=TRUE, ...) {
  s <- summary(model, ...)
  
  coefficient.names <- rownames(s$coef)
  coefficients <- s$coef[,1]
  standard.errors <- s$coef[,2]
  significance <- s$coef[,4]
  
  aic <- AIC(model)
  bic <- BIC(model)
  lik <- logLik(model)[1]
  dev <- deviance(model)
  n <- nobs(model)
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic==TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic==TRUE) {
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik==TRUE) {
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance==TRUE) {
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names=coefficient.names, 
      coef=coefficients, 
      se=standard.errors, 
      pvalues=significance, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("glm", "stats"), 
    definition = extract.glm)

extract.Relogit <- extract.glm
setMethod("extract", signature=className("Relogit", "Zelig"), 
    definition = extract.Relogit)

extract.negbin <- extract.glm
setMethod("extract", signature=className("negbin", "MASS"), 
    definition = extract.negbin)


# extension for gls objects
extract.gls <- function(model, include.aic=TRUE, include.bic=TRUE, 
    include.loglik=TRUE, include.nobs=TRUE, ...) {
  s <- summary(model, ...)
  
  coefficient.names <- rownames(s$tTable)
  coefficients <- s$tTable[,1]
  standard.errors <- s$tTable[,2]
  significance <- s$tTable[,4]
  
  lik <- s$logLik
  aic <- s$AIC
  bic <- s$BIC
  n <- nobs(model)
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic==TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic==TRUE) {
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik==TRUE) {
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names=coefficient.names, 
      coef=coefficients, 
      se=standard.errors, 
      pvalues=significance, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("gls", "nlme"), 
    definition = extract.gls)


# extension for gmm objects
extract.gmm <- function(model, include.obj.fcn=TRUE, 
    include.overidentification=FALSE, include.nobs=TRUE, ...) {
  
  s <- summary(model, ...)
  
  coefs <- s$coefficients
  names <- rownames(coefs)
  coef <- coefs[,1]
  se <- coefs[,2]
  pval <- coefs[,4]
  
  n <- model$n #number of observations
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.obj.fcn == TRUE) {
    obj.fcn <- model$objective * 10^5 #the value of the objective function
    gof <- c(gof, obj.fcn)
    gof.names <- c(gof.names, "Criterion function")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.overidentification == TRUE) {
    jtest <- s$stest$test[1]
    gof <- c(gof, jtest)
    gof.names <- c(gof.names, "J-Test")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
    coef.names=names, 
    coef=coef, 
    se=se, 
    pvalues=pval, 
    gof.names=gof.names, 
    gof=gof, 
    gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("gmm", "gmm"), 
    definition = extract.gmm)


# extension for lm objects
extract.lm <- function(model, include.rsquared=TRUE, include.adjrs=TRUE, 
    include.nobs=TRUE, ...) {
  s <- summary(model, ...)
  
  names <- rownames(s$coef)
  co <- s$coef[,1]
  se <- s$coef[,2]
  pval <- s$coef[,4]
  
  rs <- s$r.squared #extract R-squared
  adj <- s$adj.r.squared #extract adjusted R-squared
  n <- nobs(model) #extract number of observations
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.rsquared==TRUE) {
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.adjrs==TRUE) {
    gof <- c(gof, adj)
    gof.names <- c(gof.names, "Adj.\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names=names, 
      coef=co, 
      se=se, 
      pvalues=pval, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("lm", "stats"), 
    definition = extract.lm)

extract.dynlm <- extract.lm
setMethod("extract", signature=className("dynlm", "dynlm"),
    definition = extract.dynlm)


# extension for lme objects
extract.lme <- function(model, include.aic=TRUE, include.bic=TRUE, 
    include.loglik=TRUE, include.nobs=TRUE, ...) {
  s <- summary(model, ...)

  coefficient.names <- rownames(s$tTable)
  coefficients <- s$tTable[,1]
  standard.errors <- s$tTable[,2]
  significance <- s$tTable[,5]
  
  lik <- s$logLik
  aic <- s$AIC
  bic <- s$BIC
  n <- nobs(model)
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic==TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic==TRUE) {
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik==TRUE) {
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names=coefficient.names, 
      coef=coefficients, 
      se=standard.errors, 
      pvalues=significance, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("lme", "nlme"), 
    definition = extract.lme)


# extension for lmrob objects (robustbase package)
extract.lmrob <- function(model, include.nobs = TRUE, ...) {
  s <- summary(model, ...)

  names <- rownames(s$coef)
  co <- s$coef[,1]
  se <- s$coef[,2]
  pval <- s$coef[,4]
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  
  if (include.nobs == TRUE) {
    n <- length(model$residuals)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
    coef.names = names, 
    coef = co, 
    se = se,
    pvalues = pval, 
    gof.names = gof.names, 
    gof = gof, 
    gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("lmrob", "robustbase"), 
    definition = extract.lmrob)


# extension for lnam objects (sna package)
extract.lnam <- function(model, include.rsquared=TRUE, include.adjrs=TRUE, 
    include.aic=TRUE, include.bic=TRUE, include.loglik=TRUE, ...) {
  coefs <- coef(model, ...)
  coef.names <- names(coefs)
  se <- c(model$beta.se, model$rho1.se, model$rho2.se)
  p <- 2 * (1 - pnorm(abs(coefs), 0, se))
  
  rss <- sum(model$residuals^2)
  mss <- sum((model$fitted - mean(model$fitted))^2)
  rdfns <- model$df.residual + 1
  rsquared <- mss / (mss + rss)
  adj.rsquared <- 1 - (1 - mss / (mss + rss)) * model$df.total / rdfns
  lik <- model$lnlik.model
  aic <- -2 * model$lnlik.model + 2 * model$df.model
  bic <- -2 * model$lnlik.model + log(model$df.total) * model$df.model
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.rsquared==TRUE) {
    gof <- c(gof, rsquared)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.adjrs==TRUE) {
    gof <- c(gof, adj.rsquared)
    gof.names <- c(gof.names, "Adj.\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.aic==TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic==TRUE) {
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik==TRUE) {
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  
  tr <- createTexreg(
      coef.names=coef.names, 
      coef=coefs, 
      se=se, 
      pvalues=p, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("lnam", "sna"), 
    definition = extract.lnam)


# extension for lrm objects (Design or rms package); submitted by Fabrice Le Lec
extract.lrm <- function(model, include.pseudors=TRUE, include.lr=TRUE, 
    include.nobs=TRUE, ...) {
  attributes(model$coef)$names <- lapply(attributes(model$coef)$names, 
    function(x) gsub(">=", " $\\\\geq$ ", x))
  coef.names <- attributes(model$coef)$names
  coef <- model$coef
  se <- sqrt(diag(model$var))
  p <- pnorm(abs(model$coef/sqrt(diag(model$var))), 
      lower.tail = FALSE)*2
  
  pseudors <- model$stats[10] #extract pseudo R-squared
  LR <- model$stats[3] #extract LR
  n <- model$stats[1] #extract number of observations
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.pseudors==TRUE) {
    gof <- c(gof, pseudors)
    gof.names <- c(gof.names, "Pseudo R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.lr==TRUE) {
    gof <- c(gof, LR)
    gof.names <- c(gof.names, "L.R.")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names=coef.names, 
      coef=coef, 
      se=se, 
      pvalues=p, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("lrm", "rms"), 
    definition = extract.lrm)
setMethod("extract", signature=className("lrm", "Design"), 
    definition = extract.lrm)


# extension for mer (and lmerMod, glmerMod and nlmerMod) objects (lme4 package)
extract.mer <- function(model, include.pvalues=FALSE, include.aic=TRUE, 
    include.bic=TRUE, include.loglik=TRUE, include.deviance=TRUE, 
    include.nobs=TRUE, include.groups=TRUE, include.variance=TRUE, 
    mcmc.pvalues=FALSE, mcmc.size=5000, ...) {
  
  Vcov <- vcov(model, useScale = FALSE, ...)
  Vcov <- as.matrix(Vcov)
  betas <- fixef(model, ...)
  se <- sqrt(diag(Vcov))
  if (include.pvalues == FALSE) {
    #do not include p-values (default)
  } else if (mcmc.pvalues == FALSE) { #naive p-values
    zval <- betas / se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  } else if (exists("mcmcsamp")) { # mcmc-based p-values, only available up to 
    tval <- betas / se             # the stable lme4 version 0.999999-0 on CRAN;
    m <- asS4(model)               # not working with version 0.99999911-0
    ncoef <- length(betas)
    mcmc <- mcmcsamp(m, n = mcmc.size)
    hpd <- HPDinterval(mcmc)
    pval <- 2 * (1 - pt(abs(tval), nrow(m@frame) - ncoef))
  } else {
    stop(paste("MCMC sampling is not available in the currently installed", 
        "version of the lme4 package. Please use naive p-values instead or", 
        "switch off p-values completely when using texreg, or install an", 
        "older version of lme4 (like the stable CRAN release)."))
  }
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic==TRUE) {
    aic <- AIC(model)
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic==TRUE) {
    bic <- BIC(model)
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik==TRUE) {
    lik <- logLik(model)[1]
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance==TRUE) {
    dev <- deviance(model)
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs==TRUE) {
    n <- dim(model.frame(model))[1]
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.groups==TRUE) {
    grps <- sapply(model@flist, function(x) length(levels(x)))
    grp.names <- names(grps)
    grp.names <- paste("Num. groups:", grp.names)
    gof <- c(gof, grps)
    gof.names <- c(gof.names, grp.names)
    gof.decimal <- c(gof.decimal, rep(FALSE, length(grps)))
  }
  if (include.variance==TRUE) {
    vc <- VarCorr(model)
    varcomps <- c(unlist(lapply(vc, diag)),   # random intercept variances
        attr(vc, "sc")^2)                     # residual variance
    varnames <- names(varcomps)
    varnames[length(varnames)] <- "Residual"
    varnames <- paste("Variance:", varnames)
    gof <- c(gof, varcomps)
    gof.names <- c(gof.names, varnames)
    gof.decimal <- c(gof.decimal, rep(TRUE, length(varcomps)))
  }
  
  if (include.pvalues==FALSE) {
    tr <- createTexreg(
        coef.names=names(betas), 
        coef=betas, 
        se=se,
        gof.names=gof.names,
        gof=gof,
        gof.decimal=gof.decimal
    )
  } else {
    tr <- createTexreg(
        coef.names=names(betas), 
        coef=betas, 
        se=se,
        pvalues=pval,
        gof.names=gof.names,
        gof=gof,
        gof.decimal=gof.decimal
    )
  }
  return(tr)
}

setMethod("extract", signature=className("mer", "lme4"), 
    definition = extract.mer)

extract.lmerMod <- extract.mer
setMethod("extract", signature=className("lmerMod", "lme4"), 
    definition = extract.lmerMod)

extract.glmerMod <- extract.mer
setMethod("extract", signature=className("glmerMod", "lme4"), 
    definition = extract.glmerMod)

extract.nlmerMod <- extract.mer
setMethod("extract", signature=className("nlmerMod", "lme4"), 
    definition = extract.nlmerMod)


# extension for plm objects (from the plm package)
extract.plm <- function(model, include.rsquared=TRUE, include.adjrs=TRUE, 
    include.nobs=TRUE, ...) {
  s <- summary(model, ...)
  
  coefficient.names <- rownames(s$coef)
  coefficients <- s$coef[,1]
  standard.errors <- s$coef[,2]
  significance <- s$coef[,4]
  
  rs <- s$r.squared[1]
  adj <- s$r.squared[2]
  n <- length(s$resid)
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.rsquared==TRUE) {
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.adjrs==TRUE) {
    gof <- c(gof, adj)
    gof.names <- c(gof.names, "Adj.\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names=coefficient.names, 
      coef=coefficients, 
      se=standard.errors, 
      pvalues=significance, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("plm", "plm"), 
    definition = extract.plm)


# extension for pmg objects (from the plm package)
extract.pmg <- function(model, include.nobs=TRUE, ...) {
  s <- summary(model, ...)
  
  co <- s$coef
  se <- (diag(s$vcov))^(1/2) #standard errors
  t <- co / se #t-statistics
  n <- length(s$resid) #number of observations
  d <- n - length(co) #degrees of freedom
  pval <- 2 * pt(-abs(t), df=d)
  tab <- cbind(co, se, pval) #coefficient table
  names <- rownames(tab)
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names=names, 
      coef=co, 
      se=se, 
      pvalues=pval, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("pmg", "plm"), 
    definition = extract.pmg)


# extension for polr objects (MASS package)
extract.polr <- function(model, include.thresholds=FALSE, include.aic=TRUE, 
    include.bic=TRUE, include.loglik=TRUE, include.deviance=TRUE, 
    include.nobs=TRUE, ...) {
  s <- summary(model, ...)
  
  tab <- s$coefficients
  zeta.names <- names(s$zeta)
  beta <- tab[!rownames(tab) %in% zeta.names,]
  thresh <- tab[rownames(tab) %in% zeta.names,]
  
  if (sum(!rownames(tab) %in% zeta.names) == 1) {
    beta <- t(beta)
    rownames(beta) <- rownames(tab)[!rownames(tab) %in% zeta.names]
  }
  if (sum(rownames(tab) %in% zeta.names) == 1) {
    thresh <- t(thresh)
    rownames(thresh) <- rownames(tab)[rownames(tab) %in% zeta.names]
  }
  
  threshold.names <- rownames(thresh)
  threshold.coef <- thresh[,1]
  threshold.se <- thresh[,2]
  threshold.zval <- thresh[,1] / thresh[,2]
  threshold.pval <- 2 * pnorm(abs(threshold.zval), lower.tail = FALSE)
  
  beta.names <- rownames(beta)
  beta.coef <- beta[,1]
  beta.se <- beta[,2]
  beta.zval <- beta[,1] / beta[,2]
  beta.pval <- 2 * pnorm(abs(beta.zval), lower.tail = FALSE)
  
  if (include.thresholds==TRUE) {
    names <- c(beta.names, threshold.names)
    coef <- c(beta.coef, threshold.coef)
    se <- c(beta.se, threshold.se)
    pval <- c(beta.pval, threshold.pval)
  } else {
    names <- beta.names
    coef <- beta.coef
    se <- beta.se
    pval <- beta.pval
  }
  
  n <- nobs(model)
  lik <- logLik(model)[1]
  aic <- AIC(model)
  bic <- BIC(model)
  dev <- deviance(model)
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic==TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic==TRUE) {
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik==TRUE) {
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance==TRUE) {
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names=names, 
      coef=coef, 
      se=se, 
      pvalues=pval, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("polr", "MASS"), 
    definition = extract.polr)


# extension for rem.dyad objects (relevent package)
extract.rem.dyad <- function(model, include.nvertices=TRUE, 
    include.events=TRUE, include.aic=TRUE, include.aicc=TRUE, 
    include.bic=TRUE, ...) {
  
  coef <- model$coef
  coefnames <- names(coef)
  se <- diag(model$cov)^0.5
  zval <- coef/se
  pval <- 2 * (1 - pnorm(abs(zval)))
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nvertices==TRUE) {
    num.nodes <- model$n
    gof <- c(gof, num.nodes)
    gof.names <- c(gof.names, "Num.\ nodes")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.events==TRUE) {
    num.events <- model$m
    gof <- c(gof, num.events)
    gof.names <- c(gof.names, "Num.\ events")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.aic==TRUE) {
    aic <- model$AIC
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.aicc==TRUE) {
    aicc <- model$AICC
    gof <- c(gof, aicc)
    gof.names <- c(gof.names, "AICC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic==TRUE) {
    bic <- model$BIC
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  
  tr <- createTexreg(
      coef.names=coefnames, 
      coef=coef, 
      se=se, 
      pvalues=pval, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("rem.dyad", "relevent"), 
    definition = extract.rem.dyad)


# extension for rlm objects (MASS package)
extract.rlm <- function (model, include.nobs=TRUE, ...) {
  s <- summary(model, ...)
  
  names <- rownames(s$coef)
  co <- s$coef[,1]
  se <- s$coef[,2]
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    n <- nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
    coef.names = names, 
    coef = co, 
    se = se, 
    gof.names = gof.names, 
    gof = gof, 
    gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("rlm", "MASS"), 
    definition = extract.rlm)


# extension for rq objects (quantreg package)
extract.rq <- function(model, include.nobs=TRUE, include.percentile=TRUE, ...) {
  s <- summary(model, cov=TRUE, ...)
  
  co <- s$coef[,1]
  names <- rownames(s$coef)
  se <- s$coef[,2]
  pval <- s$coef[,4]
  
  n <- length(s$resid)
  tau <- s$tau
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.percentile==TRUE) {
    gof <- c(gof, tau)
    gof.names <- c(gof.names, "Percentile")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  
  tr <- createTexreg(
      coef.names=names, 
      coef=co, 
      se=se, 
      pvalues=pval, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("rq", "quantreg"), 
    definition = extract.rq)


# extension for simex objects
extract.simex <- function(model, jackknife=TRUE, include.nobs=TRUE, ...) {
  s <- summary(model, ...)
  
  if (jackknife==TRUE) {
    names <- rownames(s$coefficients$jackknife)
    co <- s$coefficients$jackknife[,1]
    se <- s$coefficients$jackknife[,2]
    pval <- s$coefficients$jackknife[,4]
  } else {
    names <- rownames(s$coefficients$asymptotic)
    co <- s$coefficients$asymptotic[,1]
    se <- s$coefficients$asymptotic[,2]
    pval <- s$coefficients$asymptotic[,4]
  }
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs==TRUE) {
    n <- length(model$model$residuals)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names=names, 
      coef=co, 
      se=se, 
      pvalues=pval, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("simex", "simex"), 
    definition = extract.simex)


# extension for stergm objects (tergm package)
extract.stergm <- function(model, beside=FALSE, include.formation=TRUE, 
    include.dissolution=TRUE, include.nvertices=TRUE, include.aic=FALSE, 
    include.bic=FALSE, include.loglik=FALSE, ...) {
  s <- summary(model, ...)
  
  if (beside==FALSE) {
    co <- numeric()
    se <- numeric()
    names <- character()
    pval <- numeric()
    if (include.formation==TRUE) {
      names <- paste("Formation:", rownames(s$formation$coefs))
      co <- s$formation$coefs[,1]
      se <- s$formation$coefs[,2]
      pval <- s$formation$coefs[,4]
    }
    if (include.dissolution==TRUE) {
      names <- c(names, paste("Dissolution:", 
          rownames(s$dissolution$coefs)))
      co <- c(co, s$dissolution$coefs[,1])
      se <- c(se, s$dissolution$coefs[,2])
      pval <- c(pval, s$dissolution$coefs[,4])
    }
    
    gof <- numeric()
    gof.names <- character()
    gof.decimal <- logical()
    if (include.nvertices==TRUE) {
      nvertices <- model$formation.fit$network$gal$n
      gof <- c(gof, nvertices)
      gof.names <- c(gof.names, "Num.\ vertices")
      gof.decimal <- c(gof.decimal, FALSE)
    }
    if (include.aic==TRUE) {
      aic.dis <- s$dissolution$aic
      aic.form <- s$formation$aic
      gof <- c(gof, aic.form, aic.dis)
      gof.names <- c(gof.names, "Formation: AIC", "Dissolution: AIC")
      gof.decimal <- c(gof.decimal, TRUE, TRUE)
    }
    if (include.bic==TRUE) {
      bic.dis <- s$dissolution$bic
      bic.form <- s$formation$bic
      gof <- c(gof, bic.form, bic.dis)
      gof.names <- c(gof.names, "Formation: BIC", "Dissolution: BIC")
      gof.decimal <- c(gof.decimal, TRUE, TRUE)
    }
    if (include.loglik==TRUE) {
      lik <- logLik(model)[1]
      gof <- c(gof, lik)
      gof.names <- c(gof.names, "Log Likelihood")
      gof.decimal <- c(gof.decimal, TRUE)
    }
    
    tr <- createTexreg(
        coef.names=names, 
        coef=co, 
        se=se, 
        pvalues=pval, 
        gof.names=gof.names, 
        gof=gof, 
        gof.decimal=gof.decimal
    )
    
    return(tr)
  } else {
    trList <- list()
    
    co <- numeric()
    se <- numeric()
    names <- character()
    pval <- numeric()
    if (include.formation==TRUE) {
      f.names <- rownames(s$formation$coefs)
      f.co <- s$formation$coefs[,1]
      f.se <- s$formation$coefs[,2]
      f.pval <- s$formation$coefs[,4]
    }
    if (include.dissolution==TRUE) {
      d.names <- rownames(s$dissolution$coefs)
      d.co <- s$dissolution$coefs[,1]
      d.se <- s$dissolution$coefs[,2]
      d.pval <- s$dissolution$coefs[,4]
    }
    
    f.gof <- numeric()
    f.gof.names <- character()
    f.gof.decimal <- logical()
    d.gof <- numeric()
    d.gof.names <- character()
    d.gof.decimal <- logical()
    if (include.nvertices==TRUE) {
      nvertices <- model$formation.fit$network$gal$n
      f.gof <- c(f.gof, nvertices)
      f.gof.names <- c(f.gof.names, "Num.\ vertices")
      f.gof.decimal <- c(f.gof.decimal, FALSE)
      d.gof <- c(d.gof, nvertices)
      d.gof.names <- c(d.gof.names, "Num.\ vertices")
      d.gof.decimal <- c(d.gof.decimal, FALSE)
    }
    if (include.aic==TRUE) {
      f.aic <- s$formation$aic
      f.gof <- c(f.gof, f.aic)
      f.gof.names <- c(f.gof.names, "AIC")
      f.gof.decimal <- c(f.gof.decimal, TRUE)
      d.aic <- s$dissolution$aic
      d.gof <- c(d.gof, d.aic)
      d.gof.names <- c(d.gof.names, "AIC")
      d.gof.decimal <- c(d.gof.decimal, TRUE)
    }
    if (include.bic==TRUE) {
      f.bic <- s$formation$bic
      f.gof <- c(f.gof, f.bic)
      f.gof.names <- c(f.gof.names, "BIC")
      f.gof.decimal <- c(f.gof.decimal, TRUE)
      d.bic <- s$dissolution$bic
      d.gof <- c(d.gof, d.bic)
      d.gof.names <- c(d.gof.names, "BIC")
      d.gof.decimal <- c(d.gof.decimal, TRUE)
    }
    if (include.loglik==TRUE) {
      lik <- logLik(model)[1]
      f.gof <- c(f.gof, lik)
      f.gof.names <- c(f.gof.names, "Log Likelihood")
      f.gof.decimal <- c(f.gof.decimal, TRUE)
      d.gof <- c(d.gof, lik)
      d.gof.names <- c(d.gof.names, "Log Likelihood")
      d.gof.decimal <- c(d.gof.decimal, TRUE)
    }
    
    if (include.formation==TRUE) {
      tr <- createTexreg(
          coef.names=f.names, 
          coef=f.co, 
          se=f.se, 
          pvalues=f.pval, 
          gof.names=f.gof.names, 
          gof=f.gof, 
          gof.decimal=f.gof.decimal
      )
      trList[[length(trList)+1]] <- tr
    }
    
    if (include.dissolution==TRUE) {
      tr <- createTexreg(
          coef.names=d.names, 
          coef=d.co, 
          se=d.se, 
          pvalues=d.pval, 
          gof.names=d.gof.names, 
          gof=d.gof, 
          gof.decimal=d.gof.decimal
      )
      trList[[length(trList)+1]] <- tr
    }
    
    return(trList)
  }
}

setMethod("extract", signature=className("stergm", "tergm"), 
    definition = extract.stergm)


# extension for svyglm objects (survey package)
extract.svyglm <- function(model, include.aic=FALSE, include.bic=FALSE, 
    include.loglik=FALSE, include.deviance=TRUE, include.dispersion=TRUE, 
    include.nobs=TRUE, ...) {
  s <- summary(model, ...)
  
  names <- rownames(coef(s))
  co <- coef(s)[,1]
  se <- coef(s)[,2]
  pval <- coef(s)[,4]
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic==TRUE) {
    aic <- AIC(model)
    if (length(aic) > 0) {
      gof <- c(gof, aic)
      gof.names <- c(gof.names, "AIC")
      gof.decimal <- c(gof.decimal, TRUE)
    } else {
      warning("AIC was not available and will be skipped!")
    }
  }
  if (include.bic==TRUE) {
    bic <- BIC(model)
    if (length(bic) > 0) {
      gof <- c(gof, bic)
      gof.names <- c(gof.names, "BIC")
      gof.decimal <- c(gof.decimal, TRUE)
    } else {
      warning("BIC was not available and will be skipped!")
    }
  }
  if (include.loglik==TRUE) {
    lik <- logLik(model)[1]
    if (length(lik) > 0) {
      gof <- c(gof, lik)
      gof.names <- c(gof.names, "Log Likelihood")
      gof.decimal <- c(gof.decimal, TRUE)
    } else {
      warning("The log likelihood was not available and will be skipped!")
    }
  }
  if (include.deviance==TRUE) {
    dev <- deviance(model)
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.dispersion==TRUE) {
    disp <- s$dispersion[1]
    gof <- c(gof, disp)
    gof.names <- c(gof.names, "Dispersion")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs==TRUE) {
    n <- nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names=names, 
      coef=co, 
      se=se, 
      pvalues=pval, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("svyglm", "survey"), 
    definition = extract.svyglm)


# extension for systemfit objects
extract.systemfit <- function(model, include.rsquared=TRUE, include.adjrs=TRUE, 
    include.nobs=TRUE, ...) {
  equationList <- list()
  for(eq in model$eq){  #go through estimated equations
    sum <- summary(eq, ...)  #extract model summary
    names <- rownames(coef(sum))
    co <- coef(sum)[,1]
    se <- coef(sum)[,2]
    pval <- coef(sum)[,4]
    
    rs <- sum$r.squared  #extract r-squared
    adj <- sum$adj.r.squared  #extract adjusted r-squared
    n <- nobs(model)  #extract number of observations
    
    gof <- numeric()
    gof.names <- character()
    gof.decimal <- logical()
    if (include.rsquared==TRUE) {
      gof <- c(gof, rs)
      gof.names <- c(gof.names, "R$^2$")
      gof.decimal <- c(gof.decimal, TRUE)
    }
    if (include.adjrs==TRUE) {
      gof <- c(gof, adj)
      gof.names <- c(gof.names, "Adj.\ R$^2$")
      gof.decimal <- c(gof.decimal, TRUE)
    }
    if (include.nobs==TRUE) {
      gof <- c(gof, n)
      gof.names <- c(gof.names, "Num.\ obs.")
      gof.decimal <- c(gof.decimal, FALSE)
    }
    
    tr <- createTexreg(
      coef.names=names, 
      coef=co, 
      se=se, 
      pvalues=pval, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
    )
    equationList[[eq$eqnNo]] <- tr
  }
  return(equationList)  #returns a list of table.content lists
}

setMethod("extract", signature=className("systemfit", "systemfit"), 
    definition = extract.systemfit)


# extension for tobit objects (AER package)
extract.tobit <- function(model, include.aic=TRUE, include.bic=TRUE, 
    include.loglik=TRUE, include.deviance=TRUE, include.nobs=FALSE, 
    include.censnobs=TRUE, include.wald=TRUE, ...) {
  s <- summary(model, ...)
  
  names <- rownames(s$coefficients)
  co <- s$coefficients[,1]
  se <- s$coefficients[,2]
  pval <- s$coefficients[,4]

  n <- nobs(model)
  censnobs <- s$n
  censnobs.names <- names(censnobs)
  aic <- AIC(model)
  bic <- BIC(model)
  lik <- logLik(model)[1]
  dev <- deviance(model)
  wald <- s$wald
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic==TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic==TRUE) {
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik==TRUE) {
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance==TRUE) {
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.censnobs==TRUE) {
    gof <- c(gof, censnobs)
    gof.names <- c(gof.names, censnobs.names)
    gof.decimal <- c(gof.decimal, rep(FALSE, length(censnobs)))
  }
  
  tr <- createTexreg(
      coef.names=names, 
      coef=co, 
      se=se, 
      pvalues=pval, 
      gof.names=gof.names, 
      gof=gof, 
      gof.decimal=gof.decimal
  )
  return(tr)
}

setMethod("extract", signature=className("tobit", "AER"), 
    definition = extract.tobit)

