# The texreg package was written by Philip Leifeld.
# Please use the forum at http://r-forge.r-project.org/projects/texreg/ 
# for bug reports, help or feature requests.


# generic extract function
setGeneric("extract", function(model, ...) standardGeneric("extract"), 
    package = "texreg")


# extension for betareg objects (betareg package)
extract.betareg <- function(model, include.precision = TRUE, 
    include.pseudors = TRUE, include.loglik = TRUE, include.nobs = TRUE, ...) {
  
  s <- summary(model, ...)
  
  coef.block <- s$coefficients$mean
  if (include.precision == TRUE) {
    phi <- s$coefficients$precision
    rownames(phi) <- paste("Precision:", rownames(phi))
    coef.block <- rbind(coef.block, phi)
  }
  names <- rownames(coef.block)
  co <- coef.block[, 1]
  se <- coef.block[, 2]
  pval <- coef.block[, 4]
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.pseudors == TRUE) {
    pseudors <- model$pseudo.r.squared
    gof <- c(gof, pseudors)
    gof.names <- c(gof.names, "Pseudo R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    lik <- model$loglik
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
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
      pvalues = pval, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("betareg", "betareg"), 
    definition = extract.betareg)


# extension for btergm objects
extract.btergm <- function(model, conf.level = 0.95, include.aic = TRUE, 
    include.bic = TRUE, include.loglik = TRUE, include.deviance = TRUE, ...) {
  #tab <- btergm.ci(model, conf.level = conf.level, print = FALSE)
  tab <- t(apply(model@bootsamp, 2, function(model) quantile(model, 
      c(((1 - conf.level) / 2), 1 - ((1 - conf.level) / 2)))))
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    aic <- AIC(model)
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    bic <- BIC(model)
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    lik <- logLik(model)
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    dev <- deviance(model)
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  
  tr <- createTexreg(
      coef.names = rownames(tab), 
      coef = coef(model), 
      ci.low = tab[, 1], 
      ci.up = tab[, 2], 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("btergm", "btergm"), 
    definition = extract.btergm)


# extension for clm objects
extract.clm <- function(model, include.thresholds = TRUE, include.aic = TRUE, 
    include.bic = TRUE, include.loglik = TRUE, include.nobs = TRUE, ...) {
  s <- summary(model, ...)
  
  tab <- s$coefficients
  thresh <- tab[rownames(tab) %in% names(s$aliased$alpha), ]
  threshold.names <- rownames(thresh)
  threshold.coef <- thresh[, 1]
  threshold.se <- thresh[, 2]
  threshold.pval <- thresh[, 4]
  beta <- tab[rownames(tab) %in% names(s$aliased$beta), ]
  beta.names <- rownames(beta)
  beta.coef <- beta[, 1]
  beta.se <- beta[, 2]
  beta.pval <- beta[, 4]
  if (include.thresholds == TRUE) {
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
  if (include.aic == TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names = names, 
      coef = coef, 
      se = se, 
      pvalues = pval, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("clm", "ordinal"), 
    definition = extract.clm)

extract.sclm <- extract.clm
setMethod("extract", signature = className("sclm", "ordinal"), 
    definition = extract.clm)


# extension for coxph objects (survival package)
extract.coxph <- function(model, include.aic = TRUE, include.rsquared = TRUE, 
    include.maxrs = TRUE, include.events = TRUE, include.nobs = TRUE, 
    include.missings = TRUE, include.zph = TRUE, ...) {
  s <- summary(model, ...)
  
  coefficient.names <- rownames(s$coef)
  coefficients <- s$coef[, 1]
  if (is.null(model$naive.var)) {
    standard.errors <- s$coef[, 3]
    significance <- s$coef[, 5]
  } else {
    standard.errors <- s$coef[, 4]
    significance <- s$coef[, 6]
  }
  
  aic <- extractAIC(model)[2]
  event <- model$nevent
  n <- model$n
  mis <- length(model$na.action)
  rs <- s$rsq[1]
  maxrs <- s$rsq[2]
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.rsquared == TRUE) {
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.maxrs == TRUE) {
    gof <- c(gof, maxrs)
    gof.names <- c(gof.names, "Max.\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.events == TRUE) {
    gof <- c(gof, event)
    gof.names <- c(gof.names, "Num.\ events")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.missings == TRUE) {
    gof <- c(gof, mis)
    gof.names <- c(gof.names, "Missings")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.zph == TRUE) {
    zph <- cox.zph(model)$table
    zph <- zph[length(zph[, 1]), length(zph[1, ])]
    gof <- c(gof, zph)
    gof.names <- c(gof.names, "PH test")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  
  tr <- createTexreg(
      coef.names = coefficient.names, 
      coef = coefficients, 
      se = standard.errors, 
      pvalues = significance, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("coxph", "survival"), 
    definition = extract.coxph)


# extension for coxph.penal objects (survival package)
extract.coxph.penal <- function(model, include.aic = TRUE, 
    include.rsquared = TRUE, include.maxrs = TRUE, include.events = TRUE, 
    include.nobs = TRUE, include.missings = TRUE, include.zph = TRUE, ...) {
  
  coefficients <- coef(model, ...)
  coefficient.names <- names(coefficients)
  standard.errors <- sqrt(diag(model$var))
  significance <- 1 - pchisq(model$coefficients^2 / diag(model$var), 1)
  
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
  if (include.aic == TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.rsquared == TRUE) {
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.maxrs == TRUE) {
    gof <- c(gof, maxrs)
    gof.names <- c(gof.names, "Max.\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.events == TRUE) {
    gof <- c(gof, event)
    gof.names <- c(gof.names, "Num.\ events")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.missings == TRUE) {
    gof <- c(gof, mis)
    gof.names <- c(gof.names, "Missings")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.zph == TRUE) {
    zph <- cox.zph(model)$table
    zph <- zph[length(zph[, 1]), length(zph[1, ])]
    gof <- c(gof, zph)
    gof.names <- c(gof.names, "PH test")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  
  tr <- createTexreg(
      coef.names = coefficient.names, 
      coef = coefficients, 
      se = standard.errors, 
      pvalues = significance, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("coxph.penal", "survival"), 
    definition = extract.coxph.penal)


# extension for clogit objects (survival package)
extract.clogit <- function(model, include.aic = TRUE, include.rsquared = TRUE, 
    include.maxrs = TRUE, include.events = TRUE, include.nobs = TRUE, 
    include.missings = TRUE, ...) {
  s <- summary(model, ...)
  
  coefficient.names <- rownames(s$coef)
  coefficients <- s$coef[, 1]
  if (is.null(model$naive.var)) {
    standard.errors <- s$coef[, 3]
    significance <- s$coef[, 5]
  } else {
    standard.errors <- s$coef[, 4]
    significance <- s$coef[, 6]
  }
  
  aic <- extractAIC(model)[2]
  event <- model$nevent
  n <- model$n
  mis <- length(model$na.action)
  rs <- s$rsq[1]
  maxrs <- s$rsq[2]
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.rsquared == TRUE) {
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.maxrs == TRUE) {
    gof <- c(gof, maxrs)
    gof.names <- c(gof.names, "Max.\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.events == TRUE) {
    gof <- c(gof, event)
    gof.names <- c(gof.names, "Num.\ events")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.missings == TRUE) {
    gof <- c(gof, mis)
    gof.names <- c(gof.names, "Missings")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names = coefficient.names, 
      coef = coefficients, 
      se = standard.errors, 
      pvalues = significance, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("clogit", "survival"), 
    definition = extract.clogit)


# extension for ergm objects
extract.ergm <- function(model, include.aic = TRUE, include.bic = TRUE, 
    include.loglik = TRUE, ...) {
  s <- summary(model, ...)
  
  coefficient.names <- rownames(s$coefs)
  coefficients <- s$coefs[, 1]
  standard.errors <- s$coefs[, 2]
  significance <- s$coefs[, 4]
  
  lik <- model$mle.lik[1]
  aic <- s$aic
  bic <- s$bic
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  
  tr <- createTexreg(
      coef.names = coefficient.names, 
      coef = coefficients, 
      se = standard.errors, 
      pvalues = significance, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("ergm", "ergm"), 
    definition = extract.ergm)


# extension for gam objects (mgcv package)
extract.gam <- function(model, include.smooth = TRUE, include.aic = TRUE, 
    include.bic = TRUE, include.loglik = TRUE, include.deviance = TRUE, 
    include.dev.expl = TRUE, include.dispersion = TRUE, include.rsquared = TRUE,
    include.gcv = TRUE, include.nobs = TRUE, include.nsmooth = TRUE, ...) {
  
  s <- summary(model, ...)
  
  coef.block <- s$p.table
  if (include.smooth == TRUE) {
    smooth <- s$s.table
    rownames(smooth) <- paste0("EDF:\ ", rownames(smooth))
    coef.block <- rbind(coef.block, smooth)
  }
  names <- rownames(coef.block)
  co <- coef.block[, 1]
  se <- coef.block[, 2]
  pval <- coef.block[, 4]
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    aic <- AIC(model)
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    bic <- BIC(model)
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    lik <- logLik(model)[1]
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log\ Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    dev <- deviance(model)
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.dev.expl == TRUE) {
    dev.expl <- s$dev.expl
    gof <- c(gof, dev.expl)
    gof.names <- c(gof.names, "Deviance\ explained")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.dispersion == TRUE) {
    dispersion <- s$dispersion
    gof <- c(gof, dispersion)
    gof.names <- c(gof.names, "Dispersion")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.rsquared == TRUE) {
    rsq <- s$r.sq
    gof <- c(gof, rsq)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.gcv == TRUE) {
    gcv <- model$gcv.ubre
    gof <- c(gof, gcv)
    gof.names <- c(gof.names, "GCV\ score")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    n <- nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.nsmooth == TRUE) {
    m <- s$m
    gof <- c(gof, m)
    gof.names <- c(gof.names, "Num.\ smooth\ terms")
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

setMethod("extract", signature = className("gam", "mgcv"), 
    definition = extract.gam)


# extension for gee objects (gee package)
extract.gee <- function(model, robust = TRUE, include.dispersion = TRUE, 
    include.nobs = TRUE, ...) {
  s <- summary(model, ...)
  
  names <- rownames(coef(s))
  co <- coef(s)[,1]
  if (robust == TRUE) {
    se <- coef(s)[, 4]
    zval <- coef(s)[, 5]
  } else {
    se <- coef(s)[, 2]
    zval <- coef(s)[, 3]
  }
  pval <- pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  
  n <- nobs(model)
  disp <- s$scale
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.dispersion == TRUE) {
    gof <- c(gof, disp)
    gof.names <- c(gof.names, "Dispersion")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
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

setMethod("extract", signature = className("gee", "gee"), 
    definition = extract.gee)


# extension for glm objects
extract.glm <- function(model, include.aic = TRUE, include.bic = TRUE, 
    include.loglik = TRUE, include.deviance = TRUE, include.nobs = TRUE, ...) {
  s <- summary(model, ...)
  
  coefficient.names <- rownames(s$coef)
  coefficients <- s$coef[, 1]
  standard.errors <- s$coef[, 2]
  significance <- s$coef[, 4]
  
  aic <- AIC(model)
  bic <- BIC(model)
  lik <- logLik(model)[1]
  dev <- deviance(model)
  n <- nobs(model)
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names = coefficient.names, 
      coef = coefficients, 
      se = standard.errors, 
      pvalues = significance, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("glm", "stats"), 
    definition = extract.glm)

extract.brglm <- extract.glm
setMethod("extract", signature = className("brglm", "brglm"), 
    definition = extract.glm)

extract.negbin <- extract.glm
setMethod("extract", signature = className("negbin", "MASS"), 
    definition = extract.negbin)


# extension for gls objects
extract.gls <- function(model, include.aic = TRUE, include.bic = TRUE, 
    include.loglik = TRUE, include.nobs = TRUE, ...) {
  s <- summary(model, ...)
  
  coefficient.names <- rownames(s$tTable)
  coefficients <- s$tTable[, 1]
  standard.errors <- s$tTable[, 2]
  significance <- s$tTable[, 4]
  
  lik <- s$logLik
  aic <- s$AIC
  bic <- s$BIC
  n <- nobs(model)
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names = coefficient.names, 
      coef = coefficients, 
      se = standard.errors, 
      pvalues = significance, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("gls", "nlme"), 
    definition = extract.gls)


# extension for gmm objects
extract.gmm <- function(model, include.obj.fcn = TRUE, 
    include.overidentification = FALSE, include.nobs = TRUE, ...) {
  
  s <- summary(model, ...)
  
  coefs <- s$coefficients
  names <- rownames(coefs)
  coef <- coefs[, 1]
  se <- coefs[, 2]
  pval <- coefs[, 4]
  
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
    coef.names = names, 
    coef = coef, 
    se = se, 
    pvalues = pval, 
    gof.names = gof.names, 
    gof = gof, 
    gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("gmm", "gmm"), 
    definition = extract.gmm)


# extension for lm objects
extract.lm <- function(model, include.rsquared = TRUE, include.adjrs = TRUE, 
    include.nobs = TRUE, ...) {
  s <- summary(model, ...)
  
  names <- rownames(s$coef)
  co <- s$coef[, 1]
  se <- s$coef[, 2]
  pval <- s$coef[, 4]
  
  rs <- s$r.squared  #extract R-squared
  adj <- s$adj.r.squared  #extract adjusted R-squared
  n <- nobs(model)  #extract number of observations
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.rsquared == TRUE) {
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.adjrs == TRUE) {
    gof <- c(gof, adj)
    gof.names <- c(gof.names, "Adj.\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
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

setMethod("extract", signature = className("lm", "stats"), 
    definition = extract.lm)

extract.dynlm <- extract.lm
setMethod("extract", signature = className("dynlm", "dynlm"),
    definition = extract.dynlm)

extract.ivreg <- extract.lm
setMethod("extract", signature = className("ivreg", "AER"),
    definition = extract.ivreg)


# extension for lme objects
extract.lme <- function(model, include.aic = TRUE, include.bic = TRUE, 
    include.loglik = TRUE, include.nobs = TRUE, ...) {
  s <- summary(model, ...)

  coefficient.names <- rownames(s$tTable)
  coefficients <- s$tTable[, 1]
  standard.errors <- s$tTable[, 2]
  significance <- s$tTable[, 5]
  
  lik <- s$logLik
  aic <- s$AIC
  bic <- s$BIC
  n <- nobs(model)
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names = coefficient.names, 
      coef = coefficients, 
      se = standard.errors, 
      pvalues = significance, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("lme", "nlme"), 
    definition = extract.lme)

extract.nlme <- extract.lme
setMethod("extract", signature = className("nlme", "nlme"), 
    definition = extract.nlme)


# extension for lmrob objects (robustbase package)
extract.lmrob <- function(model, include.nobs = TRUE, ...) {
  s <- summary(model, ...)

  names <- rownames(s$coef)
  co <- s$coef[, 1]
  se <- s$coef[, 2]
  pval <- s$coef[, 4]
  
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

setMethod("extract", signature = className("lmrob", "robustbase"), 
    definition = extract.lmrob)


# extension for lnam objects (sna package)
extract.lnam <- function(model, include.rsquared = TRUE, include.adjrs = TRUE, 
    include.aic = TRUE, include.bic = TRUE, include.loglik = TRUE, ...) {
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
  if (include.rsquared == TRUE) {
    gof <- c(gof, rsquared)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.adjrs == TRUE) {
    gof <- c(gof, adj.rsquared)
    gof.names <- c(gof.names, "Adj.\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.aic == TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  
  tr <- createTexreg(
      coef.names = coef.names, 
      coef = coefs, 
      se = se, 
      pvalues = p, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("lnam", "sna"), 
    definition = extract.lnam)


# extension for lrm objects (Design or rms package); submitted by Fabrice Le Lec
extract.lrm <- function(model, include.pseudors = TRUE, include.lr = TRUE, 
    include.nobs = TRUE, ...) {
  attributes(model$coef)$names <- lapply(attributes(model$coef)$names, 
    function(x) gsub(">=", " $\\\\geq$ ", x))
  coef.names <- attributes(model$coef)$names
  coef <- model$coef
  se <- sqrt(diag(model$var))
  p <- pnorm(abs(model$coef / sqrt(diag(model$var))), 
      lower.tail = FALSE) * 2
  
  pseudors <- model$stats[10]  # extract pseudo R-squared
  LR <- model$stats[3]  # extract LR
  n <- model$stats[1]  # extract number of observations
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.pseudors == TRUE) {
    gof <- c(gof, pseudors)
    gof.names <- c(gof.names, "Pseudo R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.lr == TRUE) {
    gof <- c(gof, LR)
    gof.names <- c(gof.names, "L.R.")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names = coef.names, 
      coef = coef, 
      se = se, 
      pvalues = p, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("lrm", "rms"), 
    definition = extract.lrm)
setMethod("extract", signature = className("lrm", "Design"), 
    definition = extract.lrm)


# extension for maBina objects (erer package)
extract.maBina <- function(model, ...) {
  
  coefficient.names <- rownames(model$out)
  coefficients <- model$out[, 1]
  standard.errors <- model$out[, 2]
  significance <- model$out[, 4]
  
  w <- extract(model$w, ...)
  gof <- w@gof
  gof.names <- w@gof.names
  gof.decimal <- w@gof.decimal
  
  tr <- createTexreg(
      coef.names = coefficient.names, 
      coef = coefficients, 
      se = standard.errors, 
      pvalues = significance, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("maBina", "erer"), 
    definition = extract.maBina)


# extension for lme4 (+ mer, lmerMod, glmerMod, nlmerMod) objects (lme4 package)
extract.lme4 <- function(model, nsim = 1000, conf.level = 0.95, 
    include.aic = TRUE, include.bic = TRUE, include.loglik = TRUE, 
    include.deviance = TRUE, include.nobs = TRUE, include.groups = TRUE, 
    include.variance = TRUE, ...) {
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    aic <- AIC(model)
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    bic <- BIC(model)
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    lik <- logLik(model)[1]
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    dev <- deviance(model)
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    n <- dim(model.frame(model))[1]
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.groups == TRUE) {
    grps <- sapply(model@flist, function(x) length(levels(x)))
    grp.names <- names(grps)
    grp.names <- paste("Num. groups:", grp.names)
    gof <- c(gof, grps)
    gof.names <- c(gof.names, grp.names)
    gof.decimal <- c(gof.decimal, rep(FALSE, length(grps)))
  }
  if (include.variance == TRUE) {
    vc <- VarCorr(model)
    varcomps <- c(unlist(lapply(vc, diag)),   # random intercept variances
        attr(vc, "sc")^2)                     # residual variance
    varnames <- names(varcomps)
    varnames[length(varnames)] <- "Residual"
    varnames <- paste("Variance:", varnames)
    if (is.na(attr(vc, "sc"))) {
      varnames <- varnames[-length(varnames)]
      varcomps <- varcomps[-length(varcomps)]
    }
    gof <- c(gof, varcomps)
    gof.names <- c(gof.names, varnames)
    gof.decimal <- c(gof.decimal, rep(TRUE, length(varcomps)))
  }
  
  betas <- fixef(model, ...)
  if (packageVersion("lme4") < 1.0) {
    if (exists("mcmcsamp")) {
      naive <- FALSE
      cat(paste0("Computing MCMC-based confidence intervals with ", nsim, 
          " iterations at a confidence level of ", conf.level, ". You may want",
          " to upgrade to a more recent version of the lme4 package.\n\n"))
      m <- asS4(model)
      mcmc <- lme4::mcmcsamp(m, n = nsim)
      hpd <- HPDinterval(mcmc, prob = conf.level)
      ci.l <- hpd$fixef[, 1]
      ci.u <- hpd$fixef[, 2]
    } else {
      naive <- TRUE
      cat(paste0("mcmcsamp function not found. Using standard errors and ", 
          "naive p values instead. Note that the p values may be inaccurate. ", 
          "You may want to upgrade to a more recent version of the lme4 ", 
          "package.\n\n"))
    }
  } else if (packageVersion("lme4") >= 1.0) {
    if ("confint.merMod" %in% methods("confint")) {
      naive <- FALSE
      cat(paste0("Computing confidence intervals with ", nsim, 
          " iterations at a confidence level of ", conf.level, ". Use ", 
          "arguments \"method = 'profile'\" or \"method = 'boot'\" ", 
          "to tweak CIs.\n\n"))
      ci <- confint(model, level = conf.level, nsim = nsim, ...)
      last <- nrow(ci)
      number <- length(betas)
      first <- last - number + 1
      ci <- ci[first:last, ]
      ci.l <- ci[, 1]
      ci.u <- ci[, 2]
    } else {
      naive <- TRUE
      cat(paste0("confint.merMod method not found. Using standard errors ", 
          "and naive p values instead. Note that the p values may be ", 
          "inaccurate. You may want to upgrade to a more recent version of ", 
          "the lme4 package.\n\n"))
    }
  }
  
  if (naive == TRUE) {
    Vcov <- vcov(model, useScale = FALSE, ...)
    Vcov <- as.matrix(Vcov)
    se <- sqrt(diag(Vcov))
    zval <- betas / se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    
    tr <- createTexreg(
        coef.names = names(betas), 
        coef = betas, 
        se = se,
        pvalues = pval,
        gof.names = gof.names,
        gof = gof,
        gof.decimal = gof.decimal
    )
  } else {
    tr <- createTexreg(
        coef.names = names(betas), 
        coef = betas, 
        ci.low = ci.l,
        ci.up = ci.u,
        gof.names = gof.names,
        gof = gof,
        gof.decimal = gof.decimal
    )
  }
  
  return(tr)
}

setMethod("extract", signature = className("lme4", "lme4"), 
    definition = extract.lme4)

extract.mer <- extract.lme4
setMethod("extract", signature = className("mer", "lme4"), 
    definition = extract.mer)

extract.lmerMod <- extract.lme4
setMethod("extract", signature = className("lmerMod", "lme4"), 
    definition = extract.lmerMod)

extract.glmerMod <- extract.lme4
setMethod("extract", signature = className("glmerMod", "lme4"), 
    definition = extract.glmerMod)

extract.nlmerMod <- extract.lme4
setMethod("extract", signature = className("nlmerMod", "lme4"), 
    definition = extract.nlmerMod)


# extension for multinom objects (nnet package)
extract.multinom <- function(model, include.pvalues = TRUE, include.aic = TRUE, 
    include.bic = TRUE, include.loglik = TRUE, include.deviance = TRUE, 
    include.nobs = TRUE, ...) {
  
  s <- summary(model, ...)
  
  names <- model$coefnames
  co <- s$coefficients
  se <- s$standard.errors
  if (include.pvalues == TRUE) {
    zval <- co / se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  } else {
    pval <- numeric(0)
  }
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    aic <- AIC(model)
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    bic <- BIC(model)
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    lik <- logLik(model)[1]
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log\ Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    dev <- deviance(model)
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    n <- length(s$fitted.values)
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

setMethod("extract", signature = className("multinom", "nnet"), 
    definition = extract.multinom)


# extension for plm objects (from the plm package)
extract.plm <- function(model, include.rsquared = TRUE, include.adjrs = TRUE, 
    include.nobs = TRUE, ...) {
  s <- summary(model, ...)
  
  coefficient.names <- rownames(s$coef)
  coefficients <- s$coef[, 1]
  standard.errors <- s$coef[, 2]
  significance <- s$coef[, 4]
  
  rs <- s$r.squared[1]
  adj <- s$r.squared[2]
  n <- length(s$resid)
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.rsquared == TRUE) {
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.adjrs == TRUE) {
    gof <- c(gof, adj)
    gof.names <- c(gof.names, "Adj.\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names = coefficient.names, 
      coef = coefficients, 
      se = standard.errors, 
      pvalues = significance, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("plm", "plm"), 
    definition = extract.plm)


# extension for pmg objects (from the plm package)
extract.pmg <- function(model, include.nobs = TRUE, ...) {
  s <- summary(model, ...)
  
  co <- s$coef
  se <- (diag(s$vcov))^(1/2)  # standard errors
  t <- co / se  # t-statistics
  n <- length(s$resid)  # number of observations
  d <- n - length(co)  # degrees of freedom
  pval <- 2 * pt(-abs(t), df = d)
  tab <- cbind(co, se, pval)  # coefficient table
  names <- rownames(tab)
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
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

setMethod("extract", signature = className("pmg", "plm"), 
    definition = extract.pmg)


# extension for polr objects (MASS package)
extract.polr <- function(model, include.thresholds = FALSE, include.aic = TRUE, 
    include.bic = TRUE, include.loglik = TRUE, include.deviance = TRUE, 
    include.nobs = TRUE, ...) {
  s <- summary(model, ...)
  
  tab <- s$coefficients
  zeta.names <- names(s$zeta)
  beta <- tab[!rownames(tab) %in% zeta.names, ]
  thresh <- tab[rownames(tab) %in% zeta.names, ]
  
  if (sum(!rownames(tab) %in% zeta.names) == 1) {
    beta <- t(beta)
    rownames(beta) <- rownames(tab)[!rownames(tab) %in% zeta.names]
  }
  if (sum(rownames(tab) %in% zeta.names) == 1) {
    thresh <- t(thresh)
    rownames(thresh) <- rownames(tab)[rownames(tab) %in% zeta.names]
  }
  
  threshold.names <- rownames(thresh)
  threshold.coef <- thresh[, 1]
  threshold.se <- thresh[, 2]
  threshold.zval <- thresh[, 1] / thresh[, 2]
  threshold.pval <- 2 * pnorm(abs(threshold.zval), lower.tail = FALSE)
  
  beta.names <- rownames(beta)
  beta.coef <- beta[, 1]
  beta.se <- beta[, 2]
  beta.zval <- beta[, 1] / beta[, 2]
  beta.pval <- 2 * pnorm(abs(beta.zval), lower.tail = FALSE)
  
  if (include.thresholds == TRUE) {
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
  if (include.aic == TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names = names, 
      coef = coef, 
      se = se, 
      pvalues = pval, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("polr", "MASS"), 
    definition = extract.polr)


# extension for rem.dyad objects (relevent package)
extract.rem.dyad <- function(model, include.nvertices = TRUE, 
    include.events = TRUE, include.aic = TRUE, include.aicc = TRUE, 
    include.bic = TRUE, ...) {
  
  coef <- model$coef
  coefnames <- names(coef)
  se <- diag(model$cov)^0.5
  zval <- coef / se
  pval <- 2 * (1 - pnorm(abs(zval)))
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nvertices == TRUE) {
    num.nodes <- model$n
    gof <- c(gof, num.nodes)
    gof.names <- c(gof.names, "Num.\ nodes")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.events == TRUE) {
    num.events <- model$m
    gof <- c(gof, num.events)
    gof.names <- c(gof.names, "Num.\ events")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.aic == TRUE) {
    aic <- model$AIC
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.aicc == TRUE) {
    aicc <- model$AICC
    gof <- c(gof, aicc)
    gof.names <- c(gof.names, "AICC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    bic <- model$BIC
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  
  tr <- createTexreg(
      coef.names = coefnames, 
      coef = coef, 
      se = se, 
      pvalues = pval, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("rem.dyad", "relevent"), 
    definition = extract.rem.dyad)


# extension for rlm objects (MASS package)
extract.rlm <- function (model, include.nobs = TRUE, ...) {
  s <- summary(model, ...)
  
  names <- rownames(s$coef)
  co <- s$coef[, 1]
  se <- s$coef[, 2]
  
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

setMethod("extract", signature = className("rlm", "MASS"), 
    definition = extract.rlm)


# extension for rq objects (quantreg package)
extract.rq <- function(model, include.nobs = TRUE, include.percentile = TRUE, 
      ...) {
  s <- summary(model, cov = TRUE, ...)
  
  co <- s$coef[, 1]
  names <- rownames(s$coef)
  se <- s$coef[, 2]
  pval <- s$coef[, 4]
  
  n <- length(s$resid)
  tau <- s$tau
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.percentile == TRUE) {
    gof <- c(gof, tau)
    gof.names <- c(gof.names, "Percentile")
    gof.decimal <- c(gof.decimal, TRUE)
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

setMethod("extract", signature = className("rq", "quantreg"), 
    definition = extract.rq)


# extension for simex objects
extract.simex <- function(model, jackknife = TRUE, include.nobs = TRUE, ...) {
  s <- summary(model, ...)
  
  if (jackknife == TRUE) {
    names <- rownames(s$coefficients$jackknife)
    co <- s$coefficients$jackknife[, 1]
    se <- s$coefficients$jackknife[, 2]
    pval <- s$coefficients$jackknife[, 4]
  } else {
    names <- rownames(s$coefficients$asymptotic)
    co <- s$coefficients$asymptotic[, 1]
    se <- s$coefficients$asymptotic[, 2]
    pval <- s$coefficients$asymptotic[, 4]
  }
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    n <- length(model$model$residuals)
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

setMethod("extract", signature = className("simex", "simex"), 
    definition = extract.simex)


# extension for stergm objects (tergm package)
extract.stergm <- function(model, beside = FALSE, include.formation = TRUE, 
    include.dissolution = TRUE, include.nvertices = TRUE, include.aic = FALSE, 
    include.bic = FALSE, include.loglik = FALSE, ...) {
  s <- summary(model, ...)
  
  if (beside == FALSE) {
    co <- numeric()
    se <- numeric()
    names <- character()
    pval <- numeric()
    if (include.formation == TRUE) {
      names <- paste("Formation:", rownames(s$formation$coefs))
      co <- s$formation$coefs[, 1]
      se <- s$formation$coefs[, 2]
      pval <- s$formation$coefs[, 4]
    }
    if (include.dissolution == TRUE) {
      names <- c(names, paste("Dissolution:", rownames(s$dissolution$coefs)))
      co <- c(co, s$dissolution$coefs[, 1])
      se <- c(se, s$dissolution$coefs[, 2])
      pval <- c(pval, s$dissolution$coefs[, 4])
    }
    
    gof <- numeric()
    gof.names <- character()
    gof.decimal <- logical()
    if (include.nvertices == TRUE) {
      nvertices <- model$formation.fit$network$gal$n
      gof <- c(gof, nvertices)
      gof.names <- c(gof.names, "Num.\ vertices")
      gof.decimal <- c(gof.decimal, FALSE)
    }
    if (include.aic == TRUE) {
      aic.dis <- s$dissolution$aic
      aic.form <- s$formation$aic
      gof <- c(gof, aic.form, aic.dis)
      gof.names <- c(gof.names, "Formation: AIC", "Dissolution: AIC")
      gof.decimal <- c(gof.decimal, TRUE, TRUE)
    }
    if (include.bic == TRUE) {
      bic.dis <- s$dissolution$bic
      bic.form <- s$formation$bic
      gof <- c(gof, bic.form, bic.dis)
      gof.names <- c(gof.names, "Formation: BIC", "Dissolution: BIC")
      gof.decimal <- c(gof.decimal, TRUE, TRUE)
    }
    if (include.loglik == TRUE) {
      lik <- logLik(model)[1]
      gof <- c(gof, lik)
      gof.names <- c(gof.names, "Log Likelihood")
      gof.decimal <- c(gof.decimal, TRUE)
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
  } else {
    trList <- list()
    
    co <- numeric()
    se <- numeric()
    names <- character()
    pval <- numeric()
    if (include.formation == TRUE) {
      f.names <- rownames(s$formation$coefs)
      f.co <- s$formation$coefs[, 1]
      f.se <- s$formation$coefs[, 2]
      f.pval <- s$formation$coefs[, 4]
    }
    if (include.dissolution == TRUE) {
      d.names <- rownames(s$dissolution$coefs)
      d.co <- s$dissolution$coefs[, 1]
      d.se <- s$dissolution$coefs[, 2]
      d.pval <- s$dissolution$coefs[, 4]
    }
    
    f.gof <- numeric()
    f.gof.names <- character()
    f.gof.decimal <- logical()
    d.gof <- numeric()
    d.gof.names <- character()
    d.gof.decimal <- logical()
    if (include.nvertices == TRUE) {
      nvertices <- model$formation.fit$network$gal$n
      f.gof <- c(f.gof, nvertices)
      f.gof.names <- c(f.gof.names, "Num.\ vertices")
      f.gof.decimal <- c(f.gof.decimal, FALSE)
      d.gof <- c(d.gof, nvertices)
      d.gof.names <- c(d.gof.names, "Num.\ vertices")
      d.gof.decimal <- c(d.gof.decimal, FALSE)
    }
    if (include.aic == TRUE) {
      f.aic <- s$formation$aic
      f.gof <- c(f.gof, f.aic)
      f.gof.names <- c(f.gof.names, "AIC")
      f.gof.decimal <- c(f.gof.decimal, TRUE)
      d.aic <- s$dissolution$aic
      d.gof <- c(d.gof, d.aic)
      d.gof.names <- c(d.gof.names, "AIC")
      d.gof.decimal <- c(d.gof.decimal, TRUE)
    }
    if (include.bic == TRUE) {
      f.bic <- s$formation$bic
      f.gof <- c(f.gof, f.bic)
      f.gof.names <- c(f.gof.names, "BIC")
      f.gof.decimal <- c(f.gof.decimal, TRUE)
      d.bic <- s$dissolution$bic
      d.gof <- c(d.gof, d.bic)
      d.gof.names <- c(d.gof.names, "BIC")
      d.gof.decimal <- c(d.gof.decimal, TRUE)
    }
    if (include.loglik == TRUE) {
      lik <- logLik(model)[1]
      f.gof <- c(f.gof, lik)
      f.gof.names <- c(f.gof.names, "Log Likelihood")
      f.gof.decimal <- c(f.gof.decimal, TRUE)
      d.gof <- c(d.gof, lik)
      d.gof.names <- c(d.gof.names, "Log Likelihood")
      d.gof.decimal <- c(d.gof.decimal, TRUE)
    }
    
    if (include.formation == TRUE) {
      tr <- createTexreg(
          coef.names = f.names, 
          coef = f.co, 
          se = f.se, 
          pvalues = f.pval, 
          gof.names = f.gof.names, 
          gof = f.gof, 
          gof.decimal = f.gof.decimal
      )
      trList[[length(trList) + 1]] <- tr
    }
    
    if (include.dissolution == TRUE) {
      tr <- createTexreg(
          coef.names = d.names, 
          coef = d.co, 
          se = d.se, 
          pvalues = d.pval, 
          gof.names = d.gof.names, 
          gof = d.gof, 
          gof.decimal = d.gof.decimal
      )
      trList[[length(trList) + 1]] <- tr
    }
    
    return(trList)
  }
}

setMethod("extract", signature = className("stergm", "tergm"), 
    definition = extract.stergm)


# extension for survreg objects (survival package)
extract.survreg <- function(model, include.aic = TRUE, include.bic = TRUE, 
    include.loglik = TRUE, include.deviance = TRUE, include.nobs = TRUE, ...) {
  
  s <- summary(model, ...)
  
  names <- rownames(s$table)
  co <- s$table[, 1]
  se <- s$table[, 2]
  pval <- s$table[, 4]
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    aic <- AIC(model)
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    bic <- BIC(model)
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    lik <- logLik(model)[1]
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log\ Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    dev <- deviance(model)
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
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
      pvalues = pval, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("survreg", "survival"), 
    definition = extract.survreg)

extract.survreg.penal <- extract.survreg
setMethod("extract", signature = className("survreg.penal", "survival"), 
    definition = extract.survreg.penal)


# extension for svyglm objects (survey package)
extract.svyglm <- function(model, include.aic = FALSE, include.bic = FALSE, 
    include.loglik = FALSE, include.deviance = TRUE, include.dispersion = TRUE, 
    include.nobs = TRUE, ...) {
  s <- summary(model, ...)
  
  names <- rownames(coef(s))
  co <- coef(s)[, 1]
  se <- coef(s)[, 2]
  pval <- coef(s)[, 4]
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    aic <- AIC(model)
    if (length(aic) > 0) {
      gof <- c(gof, aic)
      gof.names <- c(gof.names, "AIC")
      gof.decimal <- c(gof.decimal, TRUE)
    } else {
      warning("AIC was not available and will be skipped!")
    }
  }
  if (include.bic == TRUE) {
    bic <- BIC(model)
    if (length(bic) > 0) {
      gof <- c(gof, bic)
      gof.names <- c(gof.names, "BIC")
      gof.decimal <- c(gof.decimal, TRUE)
    } else {
      warning("BIC was not available and will be skipped!")
    }
  }
  if (include.loglik == TRUE) {
    lik <- logLik(model)[1]
    if (length(lik) > 0) {
      gof <- c(gof, lik)
      gof.names <- c(gof.names, "Log Likelihood")
      gof.decimal <- c(gof.decimal, TRUE)
    } else {
      warning("The log likelihood was not available and will be skipped!")
    }
  }
  if (include.deviance == TRUE) {
    dev <- deviance(model)
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.dispersion == TRUE) {
    disp <- s$dispersion[1]
    gof <- c(gof, disp)
    gof.names <- c(gof.names, "Dispersion")
    gof.decimal <- c(gof.decimal, TRUE)
  }
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
      pvalues = pval, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("svyglm", "survey"), 
    definition = extract.svyglm)


# extension for systemfit objects
extract.systemfit <- function(model, include.rsquared = TRUE, 
    include.adjrs=TRUE, include.nobs = TRUE, ...) {
  equationList <- list()
  for(eq in model$eq){  # go through estimated equations
    sum <- summary(eq, ...)  # extract model summary
    names <- rownames(coef(sum))
    co <- coef(sum)[, 1]
    se <- coef(sum)[, 2]
    pval <- coef(sum)[, 4]
    
    rs <- sum$r.squared  # extract r-squared
    adj <- sum$adj.r.squared  # extract adjusted r-squared
    n <- nobs(model)  # extract number of observations
    
    gof <- numeric()
    gof.names <- character()
    gof.decimal <- logical()
    if (include.rsquared == TRUE) {
      gof <- c(gof, rs)
      gof.names <- c(gof.names, "R$^2$")
      gof.decimal <- c(gof.decimal, TRUE)
    }
    if (include.adjrs == TRUE) {
      gof <- c(gof, adj)
      gof.names <- c(gof.names, "Adj.\ R$^2$")
      gof.decimal <- c(gof.decimal, TRUE)
    }
    if (include.nobs == TRUE) {
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
    equationList[[eq$eqnNo]] <- tr
  }
  return(equationList)  #returns a list of table.content lists
}

setMethod("extract", signature = className("systemfit", "systemfit"), 
    definition = extract.systemfit)


# extension for texreg objects (texreg package)
extract.texreg <- function(model, ...) {
  tr <- createTexreg(
    coef.names = model@coef.names,
    coef = model@coef,
    se = model@se,
    pvalues = model@pvalues,
    ci.low = model@ci.low,
    ci.up = model@ci.up,
    gof.names = model@gof.names,
    gof = model@gof,
    gof.decimal = model@gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("texreg", "texreg"), 
    definition = extract.texreg)


# extension for tobit objects (AER package)
extract.tobit <- function(model, include.aic = TRUE, include.bic = TRUE, 
    include.loglik = TRUE, include.deviance = TRUE, include.nobs = FALSE, 
    include.censnobs = TRUE, include.wald = TRUE, ...) {
  s <- summary(model, ...)
  
  names <- rownames(s$coefficients)
  co <- s$coefficients[, 1]
  se <- s$coefficients[, 2]
  pval <- s$coefficients[, 4]

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
  if (include.aic == TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.censnobs == TRUE) {
    gof <- c(gof, censnobs)
    gof.names <- c(gof.names, censnobs.names)
    gof.decimal <- c(gof.decimal, rep(FALSE, length(censnobs)))
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

setMethod("extract", signature = className("tobit", "AER"), 
    definition = extract.tobit)


# extension for weibreg objects (eha package)
extract.weibreg <- function(model, include.loglik = TRUE, include.lr = TRUE, 
    include.nobs = TRUE, include.events = TRUE, include.trisk = TRUE, ...) {
  
  coefs <- model$coefficients
  coef.names <- names(coefs)
  se <- sqrt(diag(model$var))
  pval <- 1 - pchisq((coefs / se)^2, 1)
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.loglik == TRUE) {
    lik <- model$loglik[2]
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.lr == TRUE) {
    lr <- -2 * (model$loglik[1] - model$loglik[2])
    gof <- c(gof, lr)
    gof.names <- c(gof.names, "LR test")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    n <- nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.events == TRUE) {
    ev <- model$events
    gof <- c(gof, ev)
    gof.names <- c(gof.names, "Num.\ events")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.trisk == TRUE) {
    trisk <- model$ttr
    gof <- c(gof, trisk)
    gof.names <- c(gof.names, "Total time at risk")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names = coef.names, 
      coef = coefs, 
      se = se, 
      pvalues = pval, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("weibreg", "eha"), 
    definition = extract.weibreg)

extract.phreg <- extract.weibreg
setMethod("extract", signature = className("phreg", "eha"), 
    definition = extract.phreg)

extract.aftreg <- extract.weibreg
setMethod("extract", signature = className("aftreg", "eha"), 
    definition = extract.aftreg)


# extension for zelig objects (Zelig package; tested with relogit and logit)
extract.zelig <- function(model, include.aic = TRUE, include.bic = TRUE, 
    include.loglik = TRUE, include.deviance = TRUE, include.nobs = TRUE, ...) {
  
  s <- summary(model, ...)
  
  coefficient.names <- rownames(s$coef)
  coefficients <- s$coef[, 1]
  standard.errors <- s$coef[, 2]
  significance <- s$coef[, 4]
  
  aic <- AIC(model)
  bic <- BIC(model)
  lik <- logLik(model)[1]
  dev <- s$deviance
  n <- nrow(model$data)
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names = coefficient.names, 
      coef = coefficients, 
      se = standard.errors, 
      pvalues = significance, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("zelig", "Zelig"), 
    definition = extract.zelig)


# extension for zeroinfl objects (pscl package)
extract.zeroinfl <- function(model, beside = FALSE, include.count = TRUE, 
    include.zero = TRUE, include.aic = TRUE, include.loglik = TRUE, 
    include.nobs = TRUE, ...) {
  
  s <- summary(model, ...)
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    aic <- AIC(model)
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    lik <- logLik(model)[1]
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log\ Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    n <- s$n
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  count <- s$coefficients$count
  zero <- s$coefficients$zero
  
  if (beside == FALSE) {
    if (include.count == TRUE && include.zero == TRUE) {
      rownames(count) <- paste("Count model:", rownames(count))
      rownames(zero) <- paste("Zero model:", rownames(zero))
      coef.block <- rbind(count, zero)
    } else if (include.count == TRUE) {
      coef.block <- count
    } else if (include.zero == TRUE) {
      coef.block <- zero
    } else {
      stop(paste("Either the include.count or the include.zero argument", 
          "must be TRUE."))
    }
    names <- rownames(coef.block)
    co <- coef.block[, 1]
    se <- coef.block[, 2]
    pval <- coef.block[, 4]
    
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
  } else {
    trList <- list()
    
    c.names <- rownames(count)
    c.co <- count[, 1]
    c.se <- count[, 2]
    c.pval <- count[, 4]
    z.names <- rownames(zero)
    z.co <- zero[, 1]
    z.se <- zero[, 2]
    z.pval <- zero[, 4]
    
    if (include.count == TRUE) {
      tr <- createTexreg(
          coef.names = c.names, 
          coef = c.co, 
          se = c.se, 
          pvalues = c.pval, 
          gof.names = gof.names, 
          gof = gof, 
          gof.decimal = gof.decimal
      )
      trList[[length(trList) + 1]] <- tr
    }
    
    if (include.zero == TRUE) {
      tr <- createTexreg(
          coef.names = z.names, 
          coef = z.co, 
          se = z.se, 
          pvalues = z.pval, 
          gof.names = gof.names, 
          gof = gof, 
          gof.decimal = gof.decimal
      )
      trList[[length(trList) + 1]] <- tr
    }
    if (length(trList) == 0) {
      stop(paste("Either the include.count or the include.zero argument", 
          "must be TRUE."))
    }
    return(trList)
  }
}

setMethod("extract", signature = className("zeroinfl", "pscl"), 
    definition = extract.zeroinfl)

extract.hurdle <- extract.zeroinfl
setMethod("extract", signature = className("hurdle", "pscl"), 
    definition = extract.hurdle)

