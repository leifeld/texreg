# The texreg package was written by Philip Leifeld.
# Please use the forum at http://r-forge.r-project.org/projects/texreg/ 
# for bug reports, help or feature requests.


# generic extract function
setGeneric("extract", function(model, ...) standardGeneric("extract"), 
    package="texreg")


# extension for clm objects
extract.clm <- function(model, include.thresholds=TRUE, include.aic=TRUE, 
    include.bic=TRUE, include.loglik=TRUE, include.nobs=TRUE) {

  tab <- summary(model)$coefficients
  thresh <- tab[rownames(tab) %in% names(summary(model)$aliased$alpha),]
  threshold.names <- rownames(thresh)
  threshold.coef <- thresh[,1]
  threshold.se <- thresh[,2]
  threshold.pval <- thresh[,4]
  beta <- tab[rownames(tab) %in% names(summary(model)$aliased$beta),]
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
    gof.names <- c(gof.names, "Num. obs.")
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


# extension for coxph and clogit objects (survival package)
extract.coxph <- function(model, include.aic=TRUE, include.rsquared=FALSE, 
    include.maxrs=FALSE, include.events=TRUE, include.nobs=TRUE, 
    include.missings=TRUE) {
  coefficient.names <- rownames(summary(model)$coef)
  coefficients <- summary(model)$coef[,1]
  standard.errors <- summary(model)$coef[,3]
  significance <- summary(model)$coef[,5]
  
  aic <- extractAIC(model)[2]
  event <- model$nevent
  n <- model$n
  mis <- length(model$na.action)
  rs <- summary(model)$rsq[1]
  maxrs <- summary(model)$rsq[2]
  
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
    gof.names <- c(gof.names, "Max. R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.events==TRUE) {
    gof <- c(gof, event)
    gof.names <- c(gof.names, "Num. events")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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

setMethod("extract", signature=className("coxph", "survival"), 
    definition = extract.coxph)

extract.clogit <- extract.coxph
setMethod("extract", signature=className("clogit", "survival"), 
    definition = extract.clogit)


# extension for ergm objects
extract.ergm <- function(model, include.aic=TRUE, include.bic=TRUE, 
    include.loglik=TRUE) {
  coefficient.names <- rownames(summary(model)$coefs)
  coefficients <- summary(model)$coefs[,1]
  standard.errors <- summary(model)$coefs[,2]
  significance <- summary(model)$coefs[,4]
  
  lik <- model$mle.lik[1]
  aic <- summary(model)$aic
  bic <- summary(model)$bic
  
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
    include.nobs=TRUE) {
  
  names <- rownames(coef(summary(model)))
  co <- coef(summary(model))[,1]
  if (robust==TRUE) {
    se <- coef(summary(model))[,4]
    zval <- coef(summary(model))[,5]
  } else {
    se <- coef(summary(model))[,2]
    zval <- coef(summary(model))[,3]
  }
  pval <- pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  
  n <- nobs(model)
  disp <- summary(model)$scale
  
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
    gof.names <- c(gof.names, "Num. obs.")
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
    include.loglik=TRUE, include.deviance=TRUE, include.nobs=TRUE) {
  coefficient.names <- rownames(summary(model)$coef)
  coefficients <- summary(model)$coef[,1]
  standard.errors <- summary(model)$coef[,2]
  significance <- summary(model)$coef[,4]
  
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
    gof.names <- c(gof.names, "Num. obs.")
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
    include.loglik=TRUE, include.nobs=TRUE) {
  coefficient.names <- rownames(summary(model)$tTable)
  coefficients <- summary(model)$tTable[,1]
  standard.errors <- summary(model)$tTable[,2]
  significance <- summary(model)$tTable[,4]
  
  lik <- summary(model)$logLik
  aic <- summary(model)$AIC
  bic <- summary(model)$BIC
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
    gof.names <- c(gof.names, "Num. obs.")
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


# extension for lm objects
extract.lm <- function(model, include.rsquared=TRUE, include.adjrs=TRUE, 
    include.nobs=TRUE) {
  names <- rownames(summary(model)$coef)
  co <- summary(model)$coef[,1]
  se <- summary(model)$coef[,2]
  pval <- summary(model)$coef[,4]
  
  rs <- summary(model)$r.squared #extract R-squared
  adj <- summary(model)$adj.r.squared #extract adjusted R-squared
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
    gof.names <- c(gof.names, "Adj. R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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


# extension for lme objects
extract.lme <- function(model, include.aic=TRUE, include.bic=TRUE, 
    include.loglik=TRUE, include.nobs=TRUE) {
  coefficient.names <- rownames(summary(model)$tTable)
  coefficients <- summary(model)$tTable[,1]
  standard.errors <- summary(model)$tTable[,2]
  significance <- summary(model)$tTable[,5]
  
  lik <- summary(model)$logLik
  aic <- summary(model)$AIC
  bic <- summary(model)$BIC
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
    gof.names <- c(gof.names, "Num. obs.")
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


# extension for lmerMod objects (lme4 package, version 0.99999911-0)
extract.lmerMod <- function(model, include.pvalues=FALSE, include.aic=TRUE, 
    include.bic=TRUE, include.loglik=TRUE, include.deviance=TRUE, 
    include.nobs=TRUE, include.groups=TRUE, include.variance=TRUE) {
  
  Vcov <- vcov(model, useScale = FALSE)
  Vcov <- as.matrix(Vcov)
  betas <- fixef(model)
  se <- sqrt(diag(Vcov))
  zval <- betas / se
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)

  lik <- logLik(model)[1]
  aic <- AIC(model)
  bic <- BIC(model)
  dev <- deviance(model)
  n <- dim(model.frame(model))[1]
  # n <- nobs(model)  #alternative method
  grps <- sapply(model@flist, function(x) length(levels(x)))
  grp.names <- names(grps)
  grp.names <- paste("Groups:", grp.names)
  
  vc <- VarCorr(model)
  varcomps <- c(unlist(lapply(vc, diag)),   # random intercept variances
      attr(vc, "sc")^2)                     # residual variance
  varnames <- names(varcomps)
  varnames[length(varnames)] <- "Residual"
  varnames <- gsub("\\.", "---", varnames)
  varnames <- gsub("---\\(Intercept)", "", varnames)
  varnames <- paste("Variance:", varnames)
  
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
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.groups==TRUE) {
    gof <- c(gof, grps)
    gof.names <- c(gof.names, grp.names)
    gof.decimal <- c(gof.decimal, rep(FALSE, length(grps)))
  }
  if (include.variance==TRUE) {
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

setMethod("extract", signature=className("lmerMod", "lme4"), 
    definition = extract.lmerMod)

extract.glmerMod <- extract.lmerMod
setMethod("extract", signature=className("glmerMod", "lme4"), 
    definition = extract.glmerMod)

extract.nlmerMod <- extract.lmerMod
setMethod("extract", signature=className("nlmerMod", "lme4"), 
    definition = extract.nlmerMod)


# extension for lnam objects (sna package)
extract.lnam <- function(model, include.rsquared=TRUE, include.adjrs=TRUE, 
    include.aic=TRUE, include.bic=TRUE, include.loglik=TRUE) {
  coefs <- coef(model)
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
    gof.names <- c(gof.names, "Adj. R$^2$")
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
    include.nobs=TRUE) {
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
    gof.names <- c(gof.names, "Num. obs.")
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


# extension for mer objects (lme4 package, version 0.999999-0)
extract.mer <- function(model, include.pvalues=FALSE, include.aic=TRUE, 
    include.bic=TRUE, include.loglik=TRUE, include.deviance=TRUE, 
    include.nobs=TRUE, include.groups=TRUE, include.variance=TRUE) {
  
  Vcov <- vcov(model, useScale = FALSE)
  Vcov <- as.matrix(Vcov)
  betas <- fixef(model)
  se <- sqrt(diag(Vcov))
  zval <- betas / se
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)

  lik <- logLik(model)[1]
  aic <- AIC(model)
  bic <- BIC(model)
  dev <- deviance(model)
  n <- dim(model.frame(model))[1]
  grps <- sapply(model@flist, function(x) length(levels(x)))
  grp.names <- names(grps)
  grp.names <- paste("Groups:", grp.names)
  
  vc <- VarCorr(model)
  varcomps <- c(unlist(lapply(vc, diag)),   # random intercept variances
      attr(vc, "sc")^2)                     # residual variance
  varnames <- names(varcomps)
  varnames[length(varnames)] <- "Residual"
  varnames <- gsub("\\.", "---", varnames)
  varnames <- gsub("---\\(Intercept)", "", varnames)
  varnames <- paste("Variance:", varnames)
  
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
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.groups==TRUE) {
    gof <- c(gof, grps)
    gof.names <- c(gof.names, grp.names)
    gof.decimal <- c(gof.decimal, rep(FALSE, length(grps)))
  }
  if (include.variance==TRUE) {
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


# extension for plm objects (from the plm package)
extract.plm <- function(model, include.rsquared=TRUE, include.adjrs=TRUE, 
    include.nobs=TRUE) {
  coefficient.names <- rownames(summary(model)$coef)
  coefficients <- summary(model)$coef[,1]
  standard.errors <- summary(model)$coef[,2]
  significance <- summary(model)$coef[,4]
  
  rs <- summary(model)$r.squared[1]
  adj <- summary(model)$r.squared[2]
  n <- length(summary(model)$resid)
  
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
    gof.names <- c(gof.names, "Adj. R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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
extract.pmg <- function(model, include.nobs=TRUE) {
  co <- summary(model)$coef
  se <- (diag(summary(model)$vcov))^(1/2) #standard errors
  t <- co / se #t-statistics
  n <- length(summary(model)$resid) #number of observations
  d <- n - length(co) #degrees of freedom
  pval <- 2 * pt(-abs(t), df=d)
  tab <- cbind(co, se, pval) #coefficient table
  names <- rownames(tab)
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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
extract.polr <- function(model, include.thresholds=TRUE, include.aic=TRUE, 
    include.bic=TRUE, include.loglik=TRUE, include.deviance=TRUE, 
    include.nobs=TRUE) {
  tab <- summary(model)$coefficients
  zeta.names <- names(summary(model)$zeta)
  beta <- tab[!rownames(tab) %in% zeta.names,]
  thresh <- tab[rownames(tab) %in% zeta.names,]
  
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
    gof.names <- c(gof.names, "Num. obs.")
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


# extension for rq objects (quantreg package)
extract.rq <- function(model, include.nobs=TRUE, include.percentile=TRUE) {
  co <- summary(model, cov=TRUE)$coef[,1]
  names <- rownames(summary(model, cov=TRUE)$coef)
  se <- summary(model, cov=TRUE)$coef[,2]
  pval <- summary(model, cov=TRUE)$coef[,4]
  
  n <- length(summary(model)$resid)
  tau <- summary(model)$tau
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs==TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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
extract.simex <- function(model, jackknife=TRUE, include.nobs=TRUE) {
  if (jackknife==TRUE) {
    names <- rownames(summary(model)$coefficients$jackknife)
    co <- summary(model)$coefficients$jackknife[,1]
    se <- summary(model)$coefficients$jackknife[,2]
    pval <- summary(model)$coefficients$jackknife[,4]
  } else {
    names <- rownames(summary(model)$coefficients$asymptotic)
    co <- summary(model)$coefficients$asymptotic[,1]
    se <- summary(model)$coefficients$asymptotic[,2]
    pval <- summary(model)$coefficients$asymptotic[,4]
  }
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs==TRUE) {
    n <- length(model$model$residuals)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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


# extension for svyglm objects (survey package)
extract.svyglm <- function(model, include.aic=FALSE, include.bic=FALSE, 
    include.loglik=FALSE, include.deviance=TRUE, include.dispersion=TRUE, 
    include.nobs=TRUE) {
  
  names <- rownames(coef(summary(model)))
  co <- coef(summary(model))[,1]
  se <- coef(summary(model))[,2]
  pval <- coef(summary(model))[,4]
  
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
    disp <- summary(model)$dispersion[1]
    gof <- c(gof, disp)
    gof.names <- c(gof.names, "Dispersion")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs==TRUE) {
    n <- nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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
    include.nobs=TRUE) {
  equationList <- list()
  for(eq in model$eq){  #go through estimated equations
    sum <- summary(eq)  #extract model summary
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
      gof.names <- c(gof.names, "Adj. R$^2$")
      gof.decimal <- c(gof.decimal, TRUE)
    }
    if (include.nobs==TRUE) {
      gof <- c(gof, n)
      gof.names <- c(gof.names, "Num. obs.")
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
    include.censnobs=TRUE, include.wald=TRUE) {
  
  names <- rownames(summary(model)$coefficients)
  co <- summary(model)$coefficients[,1]
  se <- summary(model)$coefficients[,2]
  pval <- summary(model)$coefficients[,4]

  n <- nobs(model)
  censnobs <- summary(model)$n
  censnobs.names <- names(censnobs)
  aic <- AIC(model)
  bic <- BIC(model)
  lik <- logLik(model)[1]
  dev <- deviance(model)
  wald <- summary(model)$wald
  
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
    gof.names <- c(gof.names, "Num. obs.")
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

