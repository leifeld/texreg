# The texreg package was written by Philip Leifeld.
# Please use the issue tracker at http://github.com/leifeld/texreg
# for bug reports, help or feature requests.


# generic extract function
setGeneric("extract", function(model, ...) standardGeneric("extract"),
    package = "texreg")

# default extract method prompts users to install the broom package
extract.broom <- function(model, ...) {
  if (!'broom' %in% row.names(installed.packages())) {
    stop("texreg does not directly support models of class ",
         class(model),
         ", but it can sometimes use the ``broom`` package to extract model information. Call texreg again after installing the ``broom`` package to see if this is possible.")
  }
  coefficients <- try(broom_coefficients(model), silent = TRUE)
  gof <- try(broom_gof(model), silent = TRUE)
  if ((class(coefficients) == 'try-error') || (class(gof) == 'try-error')) {
    stop('Neither texreg nor broom supports models of class ', class(model), '.')
  }
  tr <- createTexreg(coef.names = coefficients$term,
                     coef = coefficients$estimate,
                     se = coefficients$std.error,
                     pvalues = coefficients$p.value,
                     gof.names = gof$gof.names,
                     gof = gof$gof,
                     gof.decimal = gof$gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("ANY"),
          definition = extract.broom)

# extension for Arima objects (stats package)
extract.Arima <- function(model, include.pvalues = FALSE, include.aic = TRUE,
    include.loglik = TRUE, ...) {

  mask <- model$mask
  nam <- names(model$coef)
  co <- model$coef
  sdev <- sqrt(diag(model$var.coef))

  if (include.pvalues == TRUE) {
    t.rat <- rep(NA, length(mask))
    t.rat[mask] <- co[mask] / sdev
    pt <- 2 * pnorm(-abs(t.rat))
    setmp <- rep(NA, length(mask))
    setmp[mask] <- sdev
  } else {
    pt <- numeric()
    setmp <- sdev
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
  if (include.loglik == TRUE) {
    lik <- model$loglik
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  tr <- createTexreg(
      coef.names = nam,
      coef = co,
      se = setmp,
      pvalues = pt,
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("Arima", "stats"),
    definition = extract.Arima)


# extension for ARIMA objects (forecast package)
extract.ARIMA <- function (model, include.pvalues = FALSE, include.aic = TRUE,
    include.aicc = TRUE, include.bic = TRUE, include.loglik = TRUE, ...) {
  mask <- model$mask
  nam <- names(model$coef)
  co <- model$coef
  sdev <- sqrt(diag(model$var.coef))
  if (include.pvalues == TRUE) {
    t.rat <- rep(NA, length(mask))
    t.rat[mask] <- co[mask] / sdev
    pt <- 2 * pnorm(-abs(t.rat))
    setmp <- rep(NA, length(mask))
    setmp[mask] <- sdev
  } else {
    pt <- numeric()
    setmp <- sdev
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
  if (include.aicc == TRUE) {
    gof <- c(gof, model$aicc)
    gof.names <- c(gof.names, "AICc")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    gof <- c(gof, model$bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    lik <- model$loglik
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  tr <- createTexreg(
      coef.names = nam,
      coef = co,
      se = setmp,
      pvalues = pt,
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("ARIMA", "forecast"),
    definition = extract.ARIMA)


# extension for averaging objects (MuMIn package)
extract.averaging <- function(model, use.ci = FALSE, adjusted.se = FALSE,
    include.nobs = TRUE, ...) {

  # MuMIn >= 1.15.0 : c('coefmat.subset', 'coefmat.full')
  # MuMIn < 1.15.0 : c('coefmat', 'coefmat.full')
  ct <- summary(model)$coefmat.full
  coefs <- ct[, 1L]
  se <- ct[, if(adjusted.se) 3L else 2L]

  if (include.nobs == TRUE) {
    gof <- as.numeric(attr(model, "nobs"))
    gof.names <- "Num.\\ obs."
    gof.decimal <- FALSE
  } else {
    gof <- numeric(0L)
    gof.names <- character(0L)
    gof.decimal <- logical(0L)
  }

  if (use.ci == TRUE) {
    ci <- confint(model, full = TRUE)
    tr <- createTexreg(
      coef.names = names(coefs),
      coef = coefs,
      ci.low = ci[, 1L],
      ci.up = ci[, 2L],
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal
    )
  } else {
    tr <- createTexreg(
      coef.names = names(coefs),
      coef = coefs,
      se = se,
      pvalues = ct[, 5L],
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal
    )
  }
  return(tr)
}

setMethod("extract", signature = className("averaging", "MuMIn"),
    definition = extract.averaging)


# extension for betamfx objects (mfx package)
extract.betamfx <- function(model, include.pseudors = TRUE,
    include.loglik = TRUE, include.nobs = TRUE, ...) {
  coefnames <- rownames(model$mfxest)
  coefs <- model$mfxest[, 1]
  se <- model$mfxest[, 2]
  pval <- model$mfxest[, 4]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    gof <- c(gof, model$fit$nobs)
    gof.names <- c(gof.names, "Num.\\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, model$fit$loglik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.pseudors == TRUE) {
    gof <- c(gof, model$fit$pseudo.r.squared)
    gof.names <- c(gof.names, "Pseudo R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  tr <- createTexreg(
    coef.names = coefnames,
    coef = coefs,
    se = se,
    pvalues = pval,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("betamfx", "mfx"),
    definition = extract.betamfx)


# extension for betaor objects (mfx package)
extract.betaor <- function(model, include.pseudors = TRUE,
    include.loglik = TRUE, include.nobs = TRUE, ...) {
  coefnames <- rownames(model$oddsratio)
  coefs <- model$oddsratio[, 1]
  se <- model$oddsratio[, 2]
  pval <- model$oddsratio[, 4]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    gof <- c(gof, model$fit$nobs)
    gof.names <- c(gof.names, "Num.\\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, model$fit$loglik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.pseudors == TRUE) {
    gof <- c(gof, model$fit$pseudo.r.squared)
    gof.names <- c(gof.names, "Pseudo R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  tr <- createTexreg(
    coef.names = coefnames,
    coef = coefs,
    se = se,
    pvalues = pval,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("betaor", "mfx"),
    definition = extract.betaor)


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
extract.btergm <- function(model, level = 0.95, include.nobs = TRUE, ...) {

  tab <- btergm::confint(model, level = level)

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    gof <- c(gof, model@nobs)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }

  tr <- createTexreg(
      coef.names = rownames(tab),
      coef = tab[, 1],
      ci.low = tab[, 2],
      ci.up = tab[, 3],
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal
  )

  return(tr)
}

setMethod("extract", signature = className("btergm", "btergm"),
    definition = extract.btergm)


# extension for censReg objects (censReg package)
extract.censReg <- function(model, include.aic = TRUE, include.bic = TRUE,
    include.loglik = TRUE, include.nobs = TRUE, ...) {
  s <- summary(model, ...)

  coefs <- s$estimate[, 1]
  rn <- rownames(s$estimate)
  se <- s$estimate[, 2]
  pval <- s$estimate[, 4]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    gof <- c(gof, AIC(model)[1])
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    gof <- c(gof, BIC(model)[1])
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, logLik(model)[1])
    gof.names <- c(gof.names, "Log\ Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, s$nObs)
    gof.names <- c(gof.names, "Num.\ obs.", "Left-censored", "Uncensored",
        "Right-censored")
    gof.decimal <- c(gof.decimal, FALSE, FALSE, FALSE, FALSE)
  }

  tr <- createTexreg(
      coef.names = rn,
      coef = coefs,
      se = se,
      pvalues = pval,
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("censReg", "censReg"),
    definition = extract.censReg)


# extension for clm objects (ordinal package)
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


# extension for clmm objects (ordinal package)
extract.clmm <- function(model, include.thresholds = TRUE,
    include.loglik = TRUE, include.aic = TRUE,  include.bic = TRUE,
    include.nobs = TRUE, include.groups = TRUE, include.variance = TRUE, ...) {
  s <- summary(model, ...)

  tab <- s$coefficients
  thresh <- tab[rownames(tab) %in% names(s$alpha), , drop = FALSE]
  threshold.names <- rownames(thresh)
  threshold.coef <- thresh[, 1]
  threshold.se <- thresh[, 2]
  threshold.pval <- thresh[, 4]
  beta <- tab[rownames(tab) %in% names(s$beta), , drop = FALSE]
  beta.names <- rownames(beta)
  beta.coef <- beta[, 1]
  beta.se <- beta[, 2]
  beta.pval <- beta[, 4]

  if (include.thresholds == TRUE) {
    cfnames <- c(beta.names, threshold.names)
    coef <- c(beta.coef, threshold.coef)
    se <- c(beta.se, threshold.se)
    pval <- c(beta.pval, threshold.pval)
  } else {
    cfnames <- beta.names
    coef <- beta.coef
    se <- beta.se
    pval <- beta.pval
  }

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.loglik == TRUE) {
    lik <- logLik(model)[1]
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
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
  if (include.nobs == TRUE) {
    n <- nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.groups == TRUE) {
    grp <- s$dims$nlev.gf
    grp.names <- paste0("Groups (", names(grp), ")")
    gof <- c(gof, grp)
    gof.names <- c(gof.names, grp.names)
    gof.decimal <- c(gof.decimal, rep(FALSE, length(grp)))
  }
  if (include.variance == TRUE) {
    var.names <- character()
    var.values <- numeric()
    for (i in 1:length(s$ST)) {
      variances <- diag(s$ST[[i]] %*% t(s$ST[[i]]))
      var.names <- c(var.names, paste0("Variance: ", names(s$ST)[[i]], ": ",
          names(variances)))
      var.values <- c(var.values, variances)
    }
    gof <- c(gof, var.values)
    gof.names <- c(gof.names, var.names)
    gof.decimal <- c(gof.decimal, rep(TRUE, length(var.values)))
  }

  tr <- createTexreg(
      coef.names = cfnames,
      coef = coef,
      se = se,
      pvalues = pval,
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("clmm", "ordinal"),
    definition = extract.clmm)


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
    zph <- survival::cox.zph(model)$table
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
    zph <- survival::cox.zph(model)$table
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


# extension for coeftest objects (lmtest package)
extract.coeftest <- function(model, ...) {

  names <- rownames(model)
  co <- model[, 1]
  se <- model[, 2]
  pval <- model[, 4]

  tr <- createTexreg(
      coef.names = names,
      coef = co,
      se = se,
      pvalues = pval
  )
  return(tr)
}

setMethod("extract", signature = className("coeftest", "lmtest"),
    definition = extract.coeftest)


# extension for ergm objects
extract.ergm <- function(model, include.aic = TRUE, include.bic = TRUE,
    include.loglik = TRUE, ...) {
  s <- summary(model, ...)

  coefficient.names <- rownames(s$coefs)
  coefficients <- s$coefs[, 1]
  standard.errors <- s$coefs[, 2]
  significance <- s$coefs[, 4]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE && !is.null(s$aic)) {
    aic <- s$aic
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE && !is.null(s$bic)) {
    bic <- s$bic
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE && !is.null(model$mle.lik)) {
    lik <- model$mle.lik[1]
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


# extension for ets objects (forecast package)
extract.ets <- function (model, include.pvalues = FALSE, include.aic = TRUE,
    include.aicc = TRUE, include.bic = TRUE, include.loglik = TRUE, ...) {
  mask <- model$mask
  nam <- names(model$par)
  co <- model$par
  sdev <- rep(-Inf,length(co))
  name <- model$method
  if (include.pvalues == TRUE) {
    t.rat <- rep(NA, length(mask))
    t.rat[mask] <- co[mask] / sdev
    pt <- 2 * pnorm(-abs(t.rat))
    setmp <- rep(NA, length(mask))
    setmp[mask] <- sdev
  } else {
    pt <- numeric()
    setmp <- sdev
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
  if (include.aicc == TRUE) {
    gof <- c(gof, model$aicc)
    gof.names <- c(gof.names, "AICc")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    gof <- c(gof, model$bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    lik <- model$loglik
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  tr <- createTexreg(
      coef.names = nam,
      coef = co,
      se = setmp,
      pvalues = pt,
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal,
      model.name = name
  )
  return(tr)
}

setMethod("extract", signature = className("ets", "forecast"),
    definition = extract.ets)


# extension for ergmm objects (latentnet package)
extract.ergmm <- function(model, include.bic = TRUE, ...) {
  s <- summary(model)

  coefficient.names <- rownames(s$pmean$coef.table)
  coefficients <- s$pmean$coef.table[, 1]
  ci.low <- s$pmean$coef.table[, 2]
  ci.up <- s$pmean$coef.table[, 3]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.bic == TRUE) {
    gof <- c(gof, s$bic$overall, s$bic$Y, s$bic$Z)
    gof.names <- c(gof.names, "BIC (Overall)", "BIC (Likelihood)",
        "BIC (Latent Positions)")
    gof.decimal <- c(gof.decimal, TRUE, TRUE, TRUE)
  }

  tr <- createTexreg(
      coef.names = coefficient.names,
      coef = coefficients,
      ci.low = ci.low,
      ci.up = ci.up,
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("ergmm", "latentnet"),
    definition = extract.ergmm)


# extension for felm objects (lfe package)
extract.felm <- function(model, include.nobs = TRUE, include.rsquared = TRUE,
    include.adjrs = TRUE, include.fstatistic = FALSE, ...) {

  s <- summary(model)
  nam <- rownames(s$coefficients)
  co <- s$coefficients[, 1]
  se <- s$coefficients[, 2]
  pval <- s$coefficients[, 4]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    gof <- c(gof, s$N)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.rsquared == TRUE) {
    gof <- c(gof, s$r2, s$P.r.squared)
    gof.names <- c(gof.names, "R$^2$ (full model)", "R$^2$ (proj model)")
    gof.decimal <- c(gof.decimal, TRUE, TRUE)
  }
  if (include.adjrs == TRUE) {
    gof <- c(gof, s$r2adj, s$P.adj.r.squared)
    gof.names <- c(gof.names, "Adj.\ R$^2$ (full model)",
        "Adj.\ R$^2$ (proj model)")
    gof.decimal <- c(gof.decimal, TRUE, TRUE)
  }
  if (include.fstatistic == TRUE) {
    gof <- c(gof, s$F.fstat[1], s$F.fstat[4],
        s$P.fstat[length(s$P.fstat) - 1], s$P.fstat[1])
    gof.names <- c(gof.names, "F statistic (full model)",
        "F (full model): p-value", "F statistic (proj model)",
        "F (proj model): p-value")
    gof.decimal <- c(gof.decimal, TRUE, TRUE, TRUE, TRUE)
  }

  tr <- createTexreg(
      coef.names = nam,
      coef = co,
      se = se,
      pvalues = pval,
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("felm", "lfe"),
    definition = extract.felm)


# extension for fGARCH objects (fGarch package)
extract.fGARCH <- function(model, include.nobs = TRUE, include.aic = TRUE,
    include.loglik = TRUE, ...) {
  namesOfPars <- rownames(model@fit$matcoef)
  co <- model@fit$matcoef[, 1]
  se <- model@fit$matcoef[, 2]
  pval <- model@fit$matcoef[, 4]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    n <- length(model@data)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.aic == TRUE) {
    aic <- model@fit$ics[1]
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    lik <- model@fit$value
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  tr <- createTexreg(
      coef.names = namesOfPars,
      coef = co,
      se = se,
      pvalues = pval,
      gof.names = gof.names,
      gof.decimal = gof.decimal,
      gof = gof
  )
  return(tr)
}

setMethod("extract", signature = className("fGARCH", "fGarch"),
    definition = extract.fGARCH)


# extension for forecast objects (forecast package)
extract.forecast <- function (model, ...) {
  model <- model$model
  return(extract(model))
}

setMethod("extract", signature = className("forecast", "forecast"),
    definition = extract.forecast)


# extension for gam and bam objects (mgcv package)
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

extract.bam <- extract.gam
setMethod("extract", signature = className("bam", "mgcv"),
    definition = extract.bam)


# extension for gamlss objects (gamlss package)
extract.gamlss <- function(model, robust = FALSE, include.nobs = TRUE,
    include.nagelkerke = TRUE, include.gaic = TRUE, ...) {

  # VCOV extraction; create coefficient block
  covmat <- suppressWarnings(stats::vcov(model, type = "all", robust = robust,
      ...))
  cf <- covmat$coef  # coefficients
  namesOfPars <- names(cf)  # names of coefficients
  se <- covmat$se  # standard errors
  tvalue <- cf / se
  pvalue <-  2 * pt(-abs(tvalue), model$df.res)  # p values

  #add the parameter names to coefficients
  possiblePars <- c("$\\mu$", "$\\sigma$", "$\\nu$", "$\\tau$")
  parIndex <- 0
  parVector <- character()
  for (i in 1:length(namesOfPars)) {
    if (namesOfPars[i] == "(Intercept)") {
      parIndex <- parIndex + 1
    }
    parName <- possiblePars[parIndex]
    parVector <- c(parVector, parName)
  }
  namesOfPars <- paste(parVector, namesOfPars)

  # GOF block
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    n <- nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.nagelkerke == TRUE) {
    nk <- gamlss::Rsq(model)
    gof <- c(gof, nk)
    gof.names <- c(gof.names, "Nagelkerke R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.gaic == TRUE) {
    gaic <- gamlss::GAIC(model)
    gof <- c(gof, gaic)
    gof.names <- c(gof.names, "Generalized AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  # create and return texreg object
  tr <- createTexreg(
      coef.names = namesOfPars,
      coef = cf,
      se = se,
      pvalues = pvalue,
      gof.names = gof.names,
      gof.decimal = gof.decimal,
      gof = gof
  )
  return(tr)
}

setMethod("extract", signature = className("gamlss", "gamlss"),
    definition = extract.gamlss)


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
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)

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


# extension for geeglm objects (geepack package)
extract.geeglm <- function(model, include.scale = TRUE,
    include.correlation = TRUE, include.nobs = TRUE, ...) {
  s <- summary(model)
  names <- rownames(s$coef)
  co <- s$coef[, 1]
  se <- s$coef[, 2]
  pval <- s$coef[, 4]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()

  if (include.scale == TRUE) {
    gof = c(gof, s$geese$scale$estimate, s$geese$scale$san.se)
    gof.names = c(gof.names, "Scale parameter: gamma", "Scale parameter: SE")
    gof.decimal = c(gof.decimal, TRUE, TRUE)
  }
  if (include.correlation == TRUE) {
    gof = c(gof, s$geese$correlation$estimate, s$geese$correlation$san.se)
    gof.names = c(gof.names, "Correlation parameter: alpha",
        "Correlation parameter: SE")
    gof.decimal = c(gof.decimal, TRUE, TRUE)
  }
  if (include.nobs == TRUE) {
    n <- nrow(model.frame(model))
    nclust <- length(s$geese$clusz)
    gof = c(gof, n, nclust)
    gof.names = c(gof.names, "Num.\ obs.", "Num.\ clust.")
    gof.decimal = c(gof.decimal, FALSE, FALSE)
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

setMethod("extract", signature = className("geeglm", "geepack"),
    definition = extract.geeglm)


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


# extension for glmmadmb objects (glmmADMB package)
extract.glmmadmb <- function(model, include.variance = TRUE,
    include.dispersion = TRUE, include.zero = TRUE, include.aic = TRUE,
    include.bic = TRUE, include.loglik = TRUE, include.nobs = TRUE,
    include.groups = TRUE, ...) {

  cf <- model$b
  nam <- names(cf)
  se <- model$stdbeta
  tval <- cf / se
  pval <- 2 * pnorm(-abs(tval))

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.variance == TRUE && !is.null(model$S)) {
    vc <- nlme::VarCorr(model)
    vari <- unlist(vc)
    for (i in 1:length(vari)) {
      gof <- c(gof, vari[i])
      gof.names <- c(gof.names, paste("Variance:", names(vari)[i]))
      gof.decimal <- c(gof.decimal, TRUE)
    }
  }
  if (include.dispersion == TRUE && !is.null(model$alpha)) {
    label <- switch(model$family,
        truncnbinom = "Dispersion",
        nbinom = "Dispersion",
        gamma = "Shape",
        beta = "Dispersion",
        betabinom = "Dispersion",
        gaussian = "Residual variance",
        logistic = "Scale",
        "Dispersion"
    )
    dsp.lab <- paste0(label, ": parameter")
    sd.lab <- paste0(label, ": SD")
    disp <- model$alpha
    sd <- model$sd_alpha
    gof <- c(gof, disp, sd)
    gof.names <- c(gof.names, dsp.lab, sd.lab)
    gof.decimal <- c(gof.decimal, TRUE, TRUE)
  }
  if (include.zero == TRUE && !is.null(model$pz)) {
    zero <- model$pz
    zero.sd <- model$sd_pz
    gof <- c(gof, zero, zero.sd)
    gof.names <- c(gof.names, "Zero inflation: parameter", "Zero inflation: SD")
    gof.decimal <- c(gof.decimal, TRUE, TRUE)
  }
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
  if (include.groups == TRUE && !is.null(model$q)) {
    groups <- model$q
    for (i in 1:length(groups)) {
      gof <- c(gof, groups[i])
      gof.names <- c(gof.names, paste("Num.\ groups:", names(groups)[i]))
      gof.decimal <- c(gof.decimal, FALSE)
    }
  }

  tr <- createTexreg(
      coef.names = nam,
      coef = cf,
      se = se,
      pvalues = pval,
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("glmmadmb", "glmmADMB"),
    definition = extract.glmmadmb)


# extension for gls objects (nlme package)
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


# extension for gel objects (gmm package)
extract.gel <- function (model, include.obj.fcn = TRUE,
    include.overidentification = FALSE, include.nobs = TRUE,
    overIdentTest = c("LR", "LM", "J "), ...) {

  overIdentTest <- match.arg(overIdentTest)
  s <- summary(model, ...)
  coefs <- s$coefficients
  names <- rownames(coefs)
  coef <- coefs[, 1]
  se <- coefs[, 2]
  pval <- coefs[, 4]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.obj.fcn == TRUE) {
    obj.fcn <- model$objective * 10^5
    gof <- c(gof, obj.fcn)
    gof.names <- c(gof.names, "Criterion function")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.overidentification == TRUE) {
    w <- which(strtrim(rownames(s$stest$test), 2) == overIdentTest)
    jtest <- s$stest$test[w, ]
    gof <- c(gof, jtest)
    ntest <- rownames(s$stest$test)[w]
    gof.names <- c(gof.names, c(ntest, paste0(ntest, " p-value")))
    gof.decimal <- c(gof.decimal, TRUE, TRUE)
  }
  if (include.nobs == TRUE) {
    n <- NROW(model$gt)
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

setMethod("extract", signature = className("gel", "gmm"),
    definition = extract.gel)


# extension for gmm objects (gmm package)
extract.gmm <- function(model, include.obj.fcn = TRUE,
    include.overidentification = FALSE, include.nobs = TRUE, ...) {

  s <- summary(model, ...)

  coefs <- s$coefficients
  names <- rownames(coefs)
  coef <- coefs[, 1]
  se <- coefs[, 2]
  pval <- coefs[, 4]

  n <- model$n  # number of observations
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.obj.fcn == TRUE) {
    obj.fcn <- model$objective * 10^5  # the value of the objective function
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


# extension for H2OBinomialModel objects (h2o package)
extract.H2OBinomialModel <- function(model, standardized = FALSE,
      include.mse = TRUE, include.rsquared = TRUE, include.logloss = TRUE,
      include.meanerror = TRUE, include.auc = TRUE, include.gini = TRUE,
      include.deviance = TRUE, include.aic = TRUE, ...) {

  # extract coefficient table from model:
  coefnames <- model@model$coefficients_table$names
  if (standardized == TRUE) {
    coefs <- model@model$coefficients_table$standardized_coefficients
  } else {
    coefs <- model@model$coefficients_table$coefficients
  }

  # create empty GOF vectors and subsequently add GOF statistics from model:
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.mse == TRUE) {
    mse <- model@model$training_metrics@metrics$MSE
    gof <- c(gof, mse)
    gof.names <- c(gof.names, "MSE")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.rsquared == TRUE) {
    r2 <- model@model$training_metrics@metrics$r2
    gof <- c(gof, r2)
    gof.names <- c(gof.names, "R^2")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.logloss == TRUE) {
    logloss <- model@model$training_metrics@metrics$logloss
    gof <- c(gof, logloss)
    gof.names <- c(gof.names, "LogLoss")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.meanerror == TRUE) {
    mpce <- model@model$training_metrics@metrics$mean_per_class_error
    gof <- c(gof, mpce)
    gof.names <- c(gof.names, "Mean Per-Class Error")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.auc == TRUE) {
    auc <- model@model$training_metrics@metrics$AUC
    gof <- c(gof, auc)
    gof.names <- c(gof.names, "AUC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.gini == TRUE) {
    gini <- model@model$training_metrics@metrics$Gini
    gof <- c(gof, gini)
    gof.names <- c(gof.names, "Gini")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    nulldev <- model@model$training_metrics@metrics$null_deviance
    resdev <- model@model$training_metrics@metrics$residual_deviance
    gof <- c(gof, nulldev, resdev)
    gof.names <- c(gof.names, "Null Deviance", "Residual Deviance")
    gof.decimal <- c(gof.decimal, TRUE, TRUE)
  }
  if (include.aic == TRUE) {
    aic <- model@model$training_metrics@metrics$AIC
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  # create texreg object:
  tr <- createTexreg(
    coef.names = coefnames,
    coef = coefs,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("H2OBinomialModel", "h2o"),
    definition = extract.H2OBinomialModel)


# extension for lm objects
extract.lm <- function(model, include.rsquared = TRUE, include.adjrs = TRUE,
    include.nobs = TRUE, include.fstatistic = FALSE, include.rmse = TRUE, ...) {
  s <- summary(model, ...)

  names <- rownames(s$coefficients)
  co <- s$coefficients[, 1]
  se <- s$coefficients[, 2]
  pval <- s$coefficients[, 4]

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
  if (include.fstatistic == TRUE) {
    fstat <- s$fstatistic[[1]]
    gof <- c(gof, fstat)
    gof.names <- c(gof.names, "F statistic")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.rmse == TRUE && !is.null(s$sigma[[1]])) {
    rmse <- s$sigma[[1]]
    gof <- c(gof, rmse)
    gof.names <- c(gof.names, "RMSE")
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
    include.loglik = TRUE, include.nobs = TRUE, include.groups = TRUE,
    include.variance = FALSE, ...) {

  s <- summary(model, ...)

  coefficient.names <- rownames(s$tTable)
  coefficients <- s$tTable[, 1]
  standard.errors <- s$tTable[, 2]
  significance <- s$tTable[, 5]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    aic <- s$AIC
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    bic <- s$BIC
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    lik <- s$logLik
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
  if (include.groups == TRUE) {
    grp <- model$dims$ngrps[1:model$dims$Q]
    gof <- c(gof, grp)
    gof.names <- c(gof.names, "Num.\ groups")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.variance == TRUE ) {
    sig.all <- s$sigma
    if (!is.null(sig.all) && !is.na(sig.all)) {
      gof <- c(gof, sig.all)
      gof.names <- c(gof.names, "sigma")
      gof.decimal <- c(gof.decimal, TRUE)
    }

    vc <- nlme::VarCorr(model)
    if ("(Intercept)" %in% rownames(vc) && "StdDev" %in% colnames(vc)) {
      sig.RE <- as.numeric(vc["(Intercept)", "StdDev"])
      if (!is.null(sig.RE) && !is.na(sig.RE)) {
        gof <- c(gof, sig.RE)
        gof.names <- c(gof.names, "sigma.\ RE")
        gof.decimal <- c(gof.decimal, TRUE)
      }
    }

    cf <- coef(model$modelStruct, unconstrained = FALSE)["corStruct.Phi1"]
    rho <- unname(cf)
    if (!is.null(rho) && !is.na(rho)) {
      gof <- c(gof, rho)
      gof.names <- c(gof.names, "rho")
      gof.decimal <- c(gof.decimal, TRUE)
    }
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

extract.glmmPQL <- extract.lme
setMethod("extract", signature = className("glmmPQL", "MASS"),
    definition = extract.glmmPQL)


# extension for lme4 (+ mer, lmerMod, glmerMod, nlmerMod) objects (lme4 package)
extract.lme4 <- function(model, method = c("naive", "profile", "boot", "Wald"),
    level = 0.95, nsim = 1000, include.aic = TRUE, include.bic = TRUE,
    include.dic = FALSE, include.deviance = FALSE, include.loglik = TRUE,
    include.nobs = TRUE, include.groups = TRUE, include.variance = TRUE, ...) {

  if (packageVersion("lme4") < 1.0) {
    message("Please update to a newer 'lme4' version for full compatibility.")
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
  if (include.dic == TRUE) {  # code from the arm package, version 1.7-07
    is_REML <- lme4::isREML(model)
    llik <- logLik(model, REML = is_REML)
    dev <- deviance(lme4::refitML(model))
    n <-  lme4::getME(model, "devcomp")$dims["n"]
    Dhat <- -2 * (llik)
    pD <- dev - Dhat
    DIC <- dev + pD[[1]]
    gof <- c(gof, DIC)
    gof.names <- c(gof.names, "DIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    dev <- deviance(lme4::refitML(model))
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    lik <- logLik(model)[1]
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
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
    grp.names <- paste("Num.\ groups:", grp.names)
    gof <- c(gof, grps)
    gof.names <- c(gof.names, grp.names)
    gof.decimal <- c(gof.decimal, rep(FALSE, length(grps)))
  }
  if (include.variance == TRUE) {
    vc <- as.data.frame(lme4::VarCorr(model))
    for (i in 1:nrow(vc)) {
      if (is.na(vc[i, 2]) && is.na(vc[i, 3])) {
        gof.names <- c(gof.names, "Var: Residual")
      } else if (is.na(vc[i, 3])) {
        gof.names <- c(gof.names, paste("Var:", vc[i, 1], vc[i, 2]))
      } else {
        gof.names <- c(gof.names, paste("Cov:", vc[i, 1], vc[i, 2], vc[i, 3]))
      }
      gof <- c(gof, vc[i, 4])
      gof.decimal <- c(gof.decimal, TRUE)
    }
#    vc <- lme4::VarCorr(model)
#    varcomps <- c(unlist(lapply(vc, diag)),   # random intercept variances
#        attr(vc, "sc")^2)                     # residual variance
#    varnames <- names(varcomps)
#    varnames[length(varnames)] <- "Residual"
#    varnames <- paste("Variance:", varnames)
#    if (is.na(attr(vc, "sc"))) {
#      varnames <- varnames[-length(varnames)]
#      varcomps <- varcomps[-length(varcomps)]
#    }
#    gof <- c(gof, varcomps)
#    gof.names <- c(gof.names, varnames)
#    gof.decimal <- c(gof.decimal, rep(TRUE, length(varcomps)))
  }

  betas <- lme4::fixef(model, ...)
  if ("confint.merMod" %in% methods("confint") && method[1] != "naive") {
    ci <- tryCatch({
        ci <- confint(model, method = method[1], level = level, nsim = nsim,
        ...)
      },
      error = function(err) {
        method <- "naive"
        message(paste("Confidence intervals not available for",
            "this model. Using naive p values instead."))
      }
    )
    if (is.null(ci)) {
      method <- "naive"
    } else {
      last <- nrow(ci)
      number <- length(betas)
      first <- last - number + 1
      ci <- ci[first:last, ]
      if (class(ci) == "matrix") {
        ci.l <- ci[, 1]
        ci.u <- ci[, 2]
      } else {
        ci.l <- ci[1]
        ci.u <- ci[2]
      }
    }
  } else if (method[1] != "naive") {
    method[1] <- "naive"
    message(paste("confint.merMod method not found. Using naive p values",
        "instead."))
  }

  if (method[1] == "naive") {
    Vcov <- tryCatch({
      Vcov <- vcov(model, useScale = FALSE, ...)
    }, error = function(err) {  # Matrix package is sometimes used internally...
      stop(paste("Please load the Matrix package or update to the latest",
          "development version of lme4 and run this command again."))
    })
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


# extension for lmrob and glmrob objects (robustbase package)
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

extract.glmrob <- extract.lmrob
setMethod("extract", signature = className("glmrob", "robustbase"),
    definition = extract.glmrob)


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


# extension for logitmfx objects (mfx package)
extract.logitmfx <- function(model, include.nobs = TRUE, include.loglik = TRUE,
    include.deviance = TRUE, include.aic = TRUE, include.bic = TRUE, ...) {
  coefnames <- rownames(model$mfxest)
  coefs <- model$mfxest[, 1]
  se <- model$mfxest[, 2]
  pval <- model$mfxest[, 4]

  n <- nrow(model$fit$model)
  ll <- (model$fit$aic - (2 * length(model$fit$coefficients))) / -2

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, ll)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    gof <- c(gof, model$fit$deviance)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.aic == TRUE) {
    gof <- c(gof, model$fit$aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    bic <- (-2 * ll) + (length(model$fit$coefficients) * log(n))
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  tr <- createTexreg(
    coef.names = coefnames,
    coef = coefs,
    se = se,
    pvalues = pval,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("logitmfx", "mfx"),
    definition = extract.logitmfx)


# extension for probitmfx objects (mfx package)
extract.probitmfx <- extract.logitmfx
setMethod("extract", signature = className("probitmfx", "mfx"),
    definition = extract.probitmfx)


# extension for logitor objects (mfx package)
extract.logitor <- function(model, include.nobs = TRUE, include.loglik = TRUE,
    include.deviance = TRUE, include.aic = TRUE, include.bic = TRUE, ...) {
  coefnames <- rownames(model$oddsratio)
  coefs <- model$oddsratio[, 1]
  se <- model$oddsratio[, 2]
  pval <- model$oddsratio[, 4]

  n <- nrow(model$fit$model)
  ll <- (model$fit$aic - (2 * length(model$fit$coefficients))) / -2

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, ll)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    gof <- c(gof, model$fit$deviance)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.aic == TRUE) {
    gof <- c(gof, model$fit$aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    bic <- (-2 * ll) + (length(model$fit$coefficients) * log(n))
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  tr <- createTexreg(
    coef.names = coefnames,
    coef = coefs,
    se = se,
    pvalues = pval,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("logitor", "mfx"),
    definition = extract.logitor)


# extension for lqmm objects (lqmm package)
extract.lqmm <- function(model, include.aic = TRUE, include.bic = TRUE,
    include.loglik = TRUE, include.nobs = TRUE, include.groups = TRUE,
    include.tau = FALSE, use.ci = FALSE, beside = TRUE, ...) {

  s <- summary(model, ...)

  tau <- model$tau
  if (length(tau) == 1 && class(s$tTable) != "list") {
    tab <- list(s$tTable)  # if only one tau value, wrap in list
  } else {
    tab <- s$tTable  # multiple tau values: already wrapped in list
  }

  if (beside == TRUE) {
    trlist <- list()
    for (i in 1:length(tau)) {
      coefficient.names <- rownames(tab[[i]])
      coefficients <- tab[[i]][, 1]
      standard.errors <- tab[[i]][, 2]
      significance <- tab[[i]][, 5]
      ci.l <- tab[[i]][, 3]
      ci.u <- tab[[i]][, 4]

      gof <- numeric()
      gof.names <- character()
      gof.decimal <- logical()
      if (include.aic == TRUE) {
        gof <- c(gof, AIC(model)[i])
        gof.names <- c(gof.names, "AIC")
        gof.decimal <- c(gof.decimal, TRUE)
      }
      if (include.bic == TRUE) {
        gof <- c(gof, BIC(model)[i])
        gof.names <- c(gof.names, "BIC")
        gof.decimal <- c(gof.decimal, TRUE)
      }
      if (include.loglik == TRUE) {
        gof <- c(gof, logLik(model)[i])
        gof.names <- c(gof.names, "Log Likelihood")
        gof.decimal <- c(gof.decimal, TRUE)
      }
      if (include.nobs == TRUE) {
        n <- nobs(model)
        gof <- c(gof, n)
        gof.names <- c(gof.names, "Num.\ obs.")
        gof.decimal <- c(gof.decimal, FALSE)
      }
      if (include.groups == TRUE) {
        gof <- c(gof, model$ngroups)
        gof.names <- c(gof.names, "Num.\ groups")
        gof.decimal <- c(gof.decimal, FALSE)
      }
      if (include.tau == TRUE) {
        gof <- c(gof, tau[i])
        gof.names <- c(gof.names, "tau")
        gof.decimal <- c(gof.decimal, TRUE)
      }

      if (use.ci == FALSE) {
        tr <- createTexreg(
            coef.names = coefficient.names,
            coef = coefficients,
            se = standard.errors,
            pvalues = significance,
            gof.names = gof.names,
            gof = gof,
            gof.decimal = gof.decimal,
            model.name = as.character(tau[i])
        )
      } else {
        tr <- createTexreg(
            coef.names = coefficient.names,
            coef = coefficients,
            pvalues = significance,
            ci.low = ci.l,
            ci.up = ci.u,
            gof.names = gof.names,
            gof = gof,
            gof.decimal = gof.decimal,
            model.name = as.character(tau[i])
        )
      }
      trlist[[i]] <- tr
    }
    return(trlist)
  } else {
    coefficient.names <- character()
    coefficients <- numeric()
    standard.errors <- numeric()
    significance <- numeric()
    ci.l <- numeric()
    ci.u <- numeric()

    for (i in 1:length(tau)) {
      coefficient.names <- c(coefficient.names, paste0(rownames(tab[[i]]),
          " (", tau[i], ")"))
      coefficients <- c(coefficients, tab[[i]][, 1])
      standard.errors <- c(standard.errors, tab[[i]][, 2])
      significance <- c(significance, tab[[i]][, 5])
      ci.l <- c(ci.l, tab[[i]][, 3])
      ci.u <- c(ci.u, tab[[i]][, 4])
    }

    gof <- numeric()
    gof.names <- character()
    gof.decimal <- logical()
    if (include.aic == TRUE) {
      gof <- c(gof, AIC(model))
      gof.names <- c(gof.names, paste0("AIC (", tau, ")"))
      gof.decimal <- c(gof.decimal, rep(TRUE, length(tau)))
    }
    if (include.bic == TRUE) {
      gof <- c(gof, BIC(model))
      gof.names <- c(gof.names, paste0("BIC (", tau, ")"))
      gof.decimal <- c(gof.decimal, rep(TRUE, length(tau)))
    }
    if (include.loglik == TRUE) {
      gof <- c(gof, logLik(model))
      gof.names <- c(gof.names, paste0("Log Likelihood (", tau, ")"))
      gof.decimal <- c(gof.decimal, rep(TRUE, length(tau)))
    }
    if (include.nobs == TRUE) {
      n <- nobs(model)
      gof <- c(gof, n)
      gof.names <- c(gof.names, "Num.\ obs.")
      gof.decimal <- c(gof.decimal, FALSE)
    }
    if (include.groups == TRUE) {
      gof <- c(gof, model$ngroups)
      gof.names <- c(gof.names, "Num.\ groups")
      gof.decimal <- c(gof.decimal, FALSE)
    }

    if (use.ci == FALSE) {
      tr <- createTexreg(
          coef.names = coefficient.names,
          coef = coefficients,
          se = standard.errors,
          pvalues = significance,
          gof.names = gof.names,
          gof = gof,
          gof.decimal = gof.decimal
      )
    } else {
      tr <- createTexreg(
          coef.names = coefficient.names,
          coef = coefficients,
          pvalues = significance,
          ci.low = ci.l,
          ci.up = ci.u,
          gof.names = gof.names,
          gof = gof,
          gof.decimal = gof.decimal
      )
    }
    return(tr)
  }
}

setMethod("extract", signature = className("lqmm", "lqmm"),
    definition = extract.lqmm)



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

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    n <- model$stats[1]  # extract number of observations
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.pseudors == TRUE) {
    pseudors <- model$stats[10]  # extract pseudo R-squared
    gof <- c(gof, pseudors)
    gof.names <- c(gof.names, "Pseudo R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.lr == TRUE) {
    LR <- model$stats[3]  # extract LR
    gof <- c(gof, LR)
    gof.names <- c(gof.names, "L.R.")
    gof.decimal <- c(gof.decimal, TRUE)
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


# extension for mlogit objects (mlogit package)
extract.mlogit <- function(model, include.aic = TRUE, include.loglik = TRUE,
    include.nobs = TRUE, ...) {
  s <- summary(model, ...)

  coefs <- s$CoefTable[, 1]
  rn <- rownames(s$CoefTable)
  se <- s$CoefTable[, 2]
  pval <- s$CoefTable[, 4]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    gof <- AIC(model)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, logLik(model)[1])
    gof.names <- c(gof.names, "Log\ Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, nrow(s$residuals))
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }

  tr <- createTexreg(
      coef.names = rn,
      coef = coefs,
      se = se,
      pvalues = pval,
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("mlogit", "mlogit"),
    definition = extract.mlogit)


# extension for mnlogit objects (mnlogit package)
extract.mnlogit <- function(model, include.aic = TRUE, include.loglik = TRUE,
    include.nobs = TRUE, include.groups = TRUE, include.intercept = TRUE,
    include.iterations = FALSE, beside = FALSE, ...) {

  s <- summary(model, ...)
  coT <- s$CoefTable
  coefnames <- rownames(coT)

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    aic <- s$AIC
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    lik <- s$logLik
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log\ Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    N <- s$model.size$N
    gof <- c(gof, N)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.groups == TRUE) {
    K <- s$model.size$K
    gof <- c(gof, K)
    gof.names <- c(gof.names, "K")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.intercept == TRUE) {
    b0 <- s$model.size$intercept
    gof <- c(gof, b0)
    gof.names <- c(gof.names, "Intercept")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.iterations == TRUE) {
    iter <- s$est.stats$niters
    gradNorm <- s$est.stats$gradNorm
    diffLike <- s$est.stats$funcDiff
    gof <- c(gof, iter, gradNorm, diffLike)
    gof.names <- c(gof.names, "Iterations", "Gradient 2-norm",
        "Diff.\ Likelihood")
    gof.decimal <- c(gof.decimal, c(FALSE, TRUE, TRUE))
  }

  if (beside == FALSE) {
    co <- coT[, 1]
    se <- coT[, 2]
    pval <- coT[, 4]

    tr <- createTexreg(
        coef.names = coefnames,
        coef = co,
        se = se,
        pvalues = pval,
        gof.names = gof.names,
        gof = gof,
        gof.decimal = gof.decimal
    )
    return(tr)
  } else {
    models <- attributes(model$freq)$names[-1]
    trlist <- list()
    for (i in 1:length(models)) {
      rows <- which(grepl(paste0(models[i], "$"), coefnames))
      coeftable <- coT[rows, ]
      cn <- coefnames[rows]
      cn <- gsub(paste0(":", models[i], "$"), "", cn)
      co <- coeftable[, 1]
      se <- coeftable[, 2]
      pval <- coeftable[, 4]

      tr <- createTexreg(
          coef.names = cn,
          coef = co,
          se = se,
          pvalues = pval,
          gof.names = gof.names,
          gof = gof,
          gof.decimal = gof.decimal,
          model.name = models[i]
      )
      trlist[[i]] <- tr
    }
    return(trlist)
  }
}

setMethod("extract", signature = className("mnlogit", "mnlogit"),
    definition = extract.mnlogit)


# extension for model.selection objects (MuMIn package)
extract.model.selection <- function(model, include.loglik = TRUE,
    include.aicc = TRUE, include.delta = TRUE, include.weight = TRUE,
    include.nobs = TRUE, ...) {

  includecols <- c(loglik = include.loglik, ic = include.aicc,
      delta = include.delta, weight = include.weight)
  include <- c(includecols, nobs = include.nobs)
  decimal <- c(TRUE, TRUE, TRUE, TRUE, FALSE)[include]
  colidx <- ncol(model) - c(loglik = 3L, ic = 2L, delta = 1L, weight = 0L)
  z <- as.matrix(`[.data.frame`(model, TRUE, colidx[includecols], drop = FALSE))
  if (include.nobs) z <- cbind(z, nobs = attr(model, "nobs"))
  mode(z) <- "numeric"
  gofnames <- as.character(c(loglik = "Log Likelihood",
      ic = colnames(model)[colidx["ic"]],
      delta = "Delta", weight = "Weight",
      nobs = "Num.\\ obs.")[include])

  coeftables <- MuMIn::coefTable(model)

  ## use t-test if dfs available, otherwise z-test:
  pval <- function(ct) {
    zval <- abs(ct[, 1L] / ct[, 2L])
    2 * if (!any(is.na(ct[, 3L]))) {
      pt(zval, df = ct[, 3L], lower.tail = FALSE)
    } else {
      pnorm(zval, lower.tail = FALSE)
    }
  }

  n <- nrow(z)
  rval <- vector(length = n, mode = "list")
  for (i in 1L:n) {
    ct <- coeftables[[i]]
    rval[[i]] <- createTexreg(
      coef.names = rownames(ct),
      coef = ct[, 1L],
      se = ct[, 2L],
      pvalues = pval(ct),
      gof.names = gofnames,
      gof = z[i, ],
      gof.decimal = decimal
    )
  }
  rval
}

setMethod("extract", signature = className("model.selection", "MuMIn"),
    definition = extract.model.selection)


# extension for mtergm objects (btergm package)
extract.mtergm <- function(model, include.nobs = TRUE, include.aic = TRUE,
    include.bic = TRUE, include.loglik = TRUE, ...) {

  coefficient.names <- names(model@coef)
  coefficients <- model@coef
  standard.errors <- model@se
  significance <- model@pval

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    gof <- c(gof, model@nobs)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.aic == TRUE && !is.null(model@aic) && !is.nan(model@aic)) {
    gof <- c(gof, model@aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE && !is.null(model@bic) && !is.nan(model@bic)) {
    gof <- c(gof, model@bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE && !is.null(model@loglik) &&
      !is.nan(model@loglik)) {
    gof <- c(gof, model@loglik)
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

setMethod("extract", signature = className("mtergm", "btergm"),
    definition = extract.mtergm)


# extension for multinom objects (nnet package)
extract.multinom <- function(model, include.pvalues = TRUE, include.aic = TRUE,
    include.bic = TRUE, include.loglik = TRUE, include.deviance = TRUE,
    include.nobs = TRUE, levels = model$lev, beside = TRUE, ...) {

  s <- summary(model, ...)

  coefnames <- model$coefnames
  co <- s$coefficients
  se <- s$standard.errors

  if (class(co) != "matrix") {
    co <- t(as.matrix(co))
    se <- t(as.matrix(se))
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
    gof.names <- c(gof.names, "Log\\ Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    dev <- deviance(model)
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    n <- nrow(s$fitted.values)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }

  if (beside == TRUE) {
    trlist <- list()
    for (i in which(rownames(co) %in% levels)) {
      if (include.pvalues == TRUE) {
        zval <- co[i, ] / se[i, ]
        pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
      } else {
        pval <- numeric(0)
      }

      tr <- createTexreg(
          coef.names = coefnames,
          coef = co[i, ],
          se = se[i, ],
          pvalues = pval,
          gof.names = gof.names,
          gof = gof,
          gof.decimal = gof.decimal,
          model.name = rownames(co)[i]
      )

      trlist <- c(trlist, tr)
    }
    if (length(trlist) == 1) {
      return(trlist[[1]])
    } else {
      return(trlist)
    }
  } else {
    pval <- numeric()
    stderr <- numeric()
    coefs <- numeric()
    nm <- character()
    for (i in which(rownames(co) %in% levels)) {
      nm <- c(nm, paste0(rownames(co)[i], ": ", coefnames))
      coefs <- c(coefs, co[i, ])
      stderr <- c(stderr, se[i, ])
      if (include.pvalues == TRUE) {
        zval <- co[i, ] / se[i, ]
        pval <- c(pval, 2 * pnorm(abs(zval), lower.tail = FALSE))
      }
    }
    tr <- createTexreg(
        coef.names = nm,
        coef = coefs,
        se = stderr,
        pvalues = pval,
        gof.names = gof.names,
        gof = gof,
        gof.decimal = gof.decimal
    )
    return(tr)
  }
}

setMethod("extract", signature = className("multinom", "nnet"),
    definition = extract.multinom)


# extension for negbinirr objects (mfx package)
extract.negbinirr <- function(model, include.nobs = TRUE,
    include.loglik = TRUE, include.deviance = TRUE, include.aic = TRUE,
    include.bic = TRUE, ...) {
  coefnames <- rownames(model$irr)
  coefs <- model$irr[, 1]
  se <- model$irr[, 2]
  pval <- model$irr[, 4]

  n <- nrow(model$fit$model)
  if ("negbinirr" %in% class(model)) {
    ll <- model$fit$twologlik / 2
  } else if ("poissonirr" %in% class(model)) {
    ll <- (model$fit$aic - (2 * length(model$fit$coefficients))) / -2
  }

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, ll)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    gof <- c(gof, model$fit$deviance)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.aic == TRUE) {
    gof <- c(gof, model$fit$aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    bic <- (-2 * ll) + ((length(model$fit$coefficients) + 1) * log(n))
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  print(gof.names)
  print(gof)
  tr <- createTexreg(
    coef.names = coefnames,
    coef = coefs,
    se = se,
    pvalues = pval,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("negbinirr", "mfx"),
    definition = extract.negbinirr)

# extension for poissonirr objects (mfx package)
extract.poissonirr <- extract.negbinirr
setMethod("extract", signature = className("poissonirr", "mfx"),
    definition = extract.poissonirr)


# extension for negbinmfx objects (mfx package)
extract.negbinmfx <- function(model, include.nobs = TRUE,
    include.loglik = TRUE, include.deviance = TRUE, include.aic = TRUE,
    include.bic = TRUE, ...) {
  coefnames <- rownames(model$mfxest)
  coefs <- model$mfxest[, 1]
  se <- model$mfxest[, 2]
  pval <- model$mfxest[, 4]

  n <- nrow(model$fit$model)
  if ("negbinmfx" %in% class(model)) {
    ll <- model$fit$twologlik / 2
  } else if ("poissonmfx" %in% class(model)) {
    ll <- (model$fit$aic - (2 * length(model$fit$coefficients))) / -2
  }

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, ll)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    gof <- c(gof, model$fit$deviance)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.aic == TRUE) {
    gof <- c(gof, model$fit$aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    bic <- (-2 * ll) + ((length(model$fit$coefficients) + 1) * log(n))
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  tr <- createTexreg(
    coef.names = coefnames,
    coef = coefs,
    se = se,
    pvalues = pval,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("negbinmfx", "mfx"),
    definition = extract.negbinmfx)


# extension for poissonmfx objects (mfx package)
extract.poissonmfx <- extract.negbinmfx
setMethod("extract", signature = className("poissonmfx", "mfx"),
    definition = extract.poissonmfx)


# extension for netlogit objects (sna package)
extract.netlogit <- function(model, include.aic = TRUE, include.bic = TRUE,
    include.deviance = TRUE, include.nobs = TRUE, ...) {

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    gof <- c(gof, model$aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    gof <- c(gof, model$bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    gof <- c(gof, model$deviance, model$null.deviance)
    gof.names <- c(gof.names, "Deviance", "Null deviance")
    gof.decimal <- c(gof.decimal, TRUE, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, model$n)
    gof.names <- c(gof.names, "Num.\\ obs.")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  cf <- model$coefficients
  pv <- model$pgreqabs
  nm <- c("(Intercept)", paste0("x", 1:(length(cf) - 1)))
  if (is.null(model$dist)) {  # "classical" fit (= simple logit model)
    cvm <- chol2inv(model$qr$qr)
    se <- sqrt(diag(cvm))
  } else {  # QAP, CUG etc.
    se <- rep(NaN, length(cf))  # not perfect; results in empty brackets!
  }

  tr <- createTexreg(
      coef.names = nm,
      coef = cf,
      se = se,
      pvalues = pv,
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("netlogit", "sna"),
    definition = extract.netlogit)


# extension for ols objects (rms package)
extract.ols <- function (model, include.nobs = TRUE, include.rsquared = TRUE,
    include.adjrs = TRUE, include.fstatistic = FALSE, include.lr = TRUE, ...) {

  names <- attributes(model$coef)$names
  co <- model$coef
  se <- sqrt(diag(model$var))
  pval <- pnorm(abs(model$coef/sqrt(diag(model$var))), lower.tail = FALSE) * 2

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    n <- nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.rsquared == TRUE) {
    rs <- model$stats["R2"]
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.adjrs == TRUE) {
    adj <- 1 - (1 - model$stats["R2"]) * (n - 1) / (n - model$stats["d.f."] - 1)
    gof <- c(gof, adj)
    gof.names <- c(gof.names, "Adj.\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.fstatistic == TRUE) {
    tryCatch({
      fs <- summary.lm(model)$fstatistic[[1]]  # won't work if penalty matrix
      gof.names <- c(gof.names, "F statistic") # is given (whatever that is)
      gof.decimal <- c(gof.decimal, TRUE)
      gof <- c(gof, fs)
    }, error = {
      warning("F statistic could not be extracted.")
    }, warning = {
      warning("F statistic could not be extracted.")
    })
  }
  if (include.lr == TRUE) {
    LR <- model$stats["Model L.R."]
    gof <- c(gof, LR)
    gof.names <- c(gof.names, "L.R.")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  tr <- createTexreg(coef.names = names, coef = co, se = se, pvalues = pval,
      gof.names = gof.names, gof = gof, gof.decimal = gof.decimal)
  return(tr)
}

setMethod("extract", signature = className("ols", "rms"),
    definition = extract.ols)


# extension for pgmm objects (from the plm package)
extract.pgmm <- function(model, include.nobs = TRUE, include.sargan = TRUE,
    include.wald = TRUE, ...) {

  s <- summary(model, ...)

  coefficient.names <- rownames(s$coefficients)
  coefficients <- s$coefficients[, 1]
  standard.errors <- s$coefficients[, 2]
  significance <- s$coefficients[, 4]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    n <- attr(s, "pdim")$nT$n
    T <- attr(s, "pdim")$nT$T
    N <- attr(s, "pdim")$nT$N
    ntot <- sum(unlist(s$residuals) != 0)
    gof <- c(gof, n, T, N, ntot)
    gof.names <- c(gof.names, "n", "T", "Num.\ obs.", "Num.\ obs.\ used")
    gof.decimal <- c(gof.decimal, FALSE, FALSE, FALSE, FALSE)
  }
  if (include.sargan == TRUE) {
    sarg.stat <- s$sargan$statistic
    sarg.par <- s$sargan$parameter
    sarg.pval <- s$sargan$p.value
    gof <- c(gof, sarg.stat, sarg.par, sarg.pval)
    gof.names <- c(gof.names, "Sargan Test: chisq", "Sargan Test: df",
        "Sargan Test: p-value")
    gof.decimal <- c(gof.decimal, TRUE, TRUE, TRUE)
  }
  if (include.wald == TRUE) {
    wald.coef <- s$wald.coef$statistic[1]
    wald.pval <- s$wald.coef$p.value[1]
    wald.par <- s$wald.coef$parameter
    gof <- c(gof, wald.coef, wald.par, wald.pval)
    gof.names <- c(
        gof.names,
        "Wald Test Coefficients: chisq",
        "Wald Test Coefficients: df",
        "Wald Test Coefficients: p-value"
    )
    gof.decimal <- c(gof.decimal, TRUE, FALSE, TRUE)
    if (!is.null(s$wald.td)) {
      td.coef <- s$wald.td$statistic[1]
      td.pval <- s$wald.td$p.value[1]
      td.par <- s$wald.td$parameter
      gof <- c(gof, td.coef, td.par, td.pval)
      gof.names <- c(
          gof.names,
          "Wald Test Time Dummies: chisq",
          "Wald Test Time Dummies: df",
          "Wald Test Time Dummies: p-value"
      )
      gof.decimal <- c(gof.decimal, TRUE, FALSE, TRUE)
    }
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

setMethod("extract", signature = className("pgmm", "plm"),
    definition = extract.pgmm)


# extension for plm objects (from the plm package)
extract.plm <- function(model, include.rsquared = TRUE, include.adjrs = TRUE,
    include.nobs = TRUE, include.variance = TRUE, ...) {
  s <- summary(model, ...)

  coefficient.names <- rownames(coef(s))
  coefficients <- coef(s)[, 1]
  standard.errors <- coef(s)[, 2]
  significance <- coef(s)[, 4]

  rs <- s$r.squared[1]
  adj <- s$r.squared[2]
  n <- length(model$residuals)

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.variance == TRUE) {
    if (model$args$model == "random") {
      se <- sqrt(unlist(plm::ercomp(model)$sigma2))
      gof <- c(gof, se)
      gof.names <- c(gof.names, paste0("s_", names(se)))
      gof.decimal <- c(gof.decimal, rep(TRUE, length(se)))
    }
  }
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


# extension for sarlm objects (spdep package)
extract.sarlm <- function(model, include.nobs = TRUE, include.loglik = TRUE,
    include.aic = TRUE, include.lr = TRUE, include.wald = TRUE, ...) {
  s <- summary(model, ...)

  names <- rownames(s$Coef)
  cf <- s$Coef[, 1]
  se <- s$Coef[, 2]
  p <- s$Coef[, ncol(s$Coef)]

  if (model$type != "error") {  # include coefficient for autocorrelation term
    rho <- model$rho
    cf <- c(cf, rho)
    names <- c(names, "$\\rho$")
    if (!is.null(model$rho.se)) {
      if (!is.null(model$adj.se)) {
        rho.se <- sqrt((model$rho.se^2) * model$adj.se)
      } else {
        rho.se <- model$rho.se
      }
      rho.pval <- 2 * (1 - pnorm(abs(rho / rho.se)))
      se <- c(se, rho.se)
      p <- c(p, rho.pval)
    } else {
      se <- c(se, NA)
      p <- c(p, NA)
    }
  }

  if (!is.null(model$lambda)) {
    cf <-c(cf, model$lambda)
    names <- c(names, "$\\lambda$")
    if (!is.null(model$lambda.se)) {
      if (!is.null(model$adj.se)) {
        lambda.se <- sqrt((model$lambda.se^2) * model$adj.se)
      } else {
        lambda.se <- model$lambda.se
      }
      lambda.pval <- 2 * (1 - pnorm(abs(model$lambda / lambda.se)))
      se <- c(se, lambda.se)
      p <- c(p, lambda.pval)
    } else {
      se <- c(se, NA)
      p <- c(p, NA)
    }
  }

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()

  if (include.nobs == TRUE) {
    n <- length(s$fitted.values)
    param <- s$parameters
    gof <- c(gof, n, param)
    gof.names <- c(gof.names, "Num.\ obs.", "Parameters")
    gof.decimal <- c(gof.decimal, FALSE, FALSE)
  }
  if (include.loglik == TRUE) {
    ll <- s$LL
    gof <- c(gof, ll)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.aic == TRUE) {
    aic <- AIC(model)
    aiclm <- s$AIC_lm.model
    gof <- c(gof, aiclm, aic)
    gof.names <- c(gof.names, "AIC (Linear model)", "AIC (Spatial model)")
    gof.decimal <- c(gof.decimal, TRUE, TRUE)
  }
  if (include.lr == TRUE && !is.null(s$LR1)) {
    gof <- c(gof, s$LR1$statistic[[1]], s$LR1$p.value[[1]])
    gof.names <- c(gof.names, "LR test: statistic", "LR test: p-value")
    gof.decimal <- c(gof.decimal, TRUE, TRUE)
  }
  if (include.wald == TRUE && !is.null(model$Wald1)) {
    waldstat <- model$Wald1$statistic
    waldp <- model$Wald1$p.value
    gof <- c(gof, waldstat, waldp)
    gof.names <- c(gof.names, "Wald test: statistic", "Wald test: p-value")
    gof.decimal <- c(gof.decimal, TRUE, TRUE)
  }

  tr <- createTexreg(
      coef.names = names,
      coef = cf,
      se = se,
      pvalues = p,
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("sarlm", "spdep"),
    definition = extract.sarlm)


# extension for selection objects (sampleSelection package)
extract.selection <- function(model, prefix = TRUE, include.selection = TRUE,
    include.outcome = TRUE, include.errors = TRUE, include.aic = TRUE,
    include.bic = TRUE, include.loglik = TRUE, include.rsquared = TRUE,
    include.adjrs = TRUE, include.nobs = TRUE, ...) {

  # extract coefficients etc.
  s <- summary(model, ...)
  coefs <- coef(model)
  coef.tab <- coef(s)
  rn <- names(coefs)
  se <- coef.tab[, 2]
  p <- coef.tab[, 4]

  # add prefixes to labels of selection and outcome components
  indices.selection <- s$param$index$betaS
  if (model$tobitType == 5) {
	  indices.outcome <- s$param$index$outcome
  } else if(model$tobitType == 2) {
		indices.outcome <- s$param$index$betaO
  }
  indices.errorterms <- s$param$index$errTerms
  if (prefix == TRUE) {
    rn[indices.selection] <- paste("S:", rn[indices.selection])
    rn[indices.outcome] <- paste("O:", rn[indices.outcome])
  }

  # retain only those components that are set up in the arguments
  include <- numeric()
  if (include.selection == TRUE) {
    include <- c(include, indices.selection)
  }
  if (include.outcome == TRUE) {
    include <- c(include, indices.outcome)
  }
  if (include.errors == TRUE) {
    include <- c(include, indices.errorterms)
  }
  coefs <- coefs[include]
  rn <- rn[include]
  se <- se[include]
  p <- p[include]

  # GOF block
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE && "loglik" %in% names(s)) {
    gof <- c(gof, AIC(model))
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE && "loglik" %in% names(s)) {
    gof <- c(gof, BIC(model))
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE && "loglik" %in% names(s)) {
    gof <- c(gof, logLik(model)[1])
    gof.names <- c(gof.names, "Log\ Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.rsquared == TRUE && "rSquared" %in% names(s)) {
    gof <- c(gof, s$rSquared$R2)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.adjrs == TRUE && "rSquared" %in% names(s)) {
    gof <- c(gof, s$rSquared$R2adj)
    gof.names <- c(gof.names, "Adj.\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    if(model$tobitType == 5) {
	    gof <- c(gof, s$param$nObs, s$param$N1, s$param$N2)
	  } else if(model$tobitType == 2) {
	    gof <- c(gof, s$param$nObs, s$param$N0, s$param$N1)
	  }
    gof.names <- c(gof.names, "Num.\ obs.", "Censored", "Observed")
    gof.decimal <- c(gof.decimal, FALSE, FALSE, FALSE)
  }

  # create and return texreg object
  tr <- createTexreg(
    coef.names = rn,
    coef = coefs,
    se = se,
    pvalues = p,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
}

setMethod("extract", signature = className("selection", "sampleSelection"),
    definition = extract.selection)


# extension for sienaFit objects (RSiena package)
extract.sienaFit <- function(model, include.iterations = TRUE, ...) {

  s <- summary(model, ...)

  coefs <- c(model$rate, model$theta)

  theta.names <- s$effects$effectName
  if (length(theta.names) == length(coefs)) {
    coef.names <- theta.names
  } else {
    if (!is.null(model$rate)) {
      rate <- model$rate
    } else {
      rate <- which(model$effects$type == "rate")
    }
    rate.names <- paste("Rate parameter period", 1:length(rate))
    coef.names <- c(rate.names, theta.names)
  }

  se <- c(model$vrate, sqrt(diag(model$covtheta)))

  pval <- 2 * pnorm(-abs(coefs / se))

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.iterations == TRUE) {
    n <- s$n
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Iterations")
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

setMethod("extract", signature = className("sienaFit", "RSiena"),
    definition = extract.sienaFit)


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
          gof.decimal = f.gof.decimal,
          model.name = "Formation"
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
          gof.decimal = d.gof.decimal,
          model.name = "Dissolution"
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
  pval <- s$table[, ncol(s$table)]

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
    if (!is.null(bic) && !is.na(bic)) {
      gof <- c(gof, bic)
      gof.names <- c(gof.names, "BIC")
      gof.decimal <- c(gof.decimal, TRUE)
    }
  }
  if (include.loglik == TRUE) {
    lik <- logLik(model)[1]
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    dev <- deviance(model)
    if (!is.null(dev)) {
      gof <- c(gof, dev)
      gof.names <- c(gof.names, "Deviance")
      gof.decimal <- c(gof.decimal, TRUE)
    }
  }
  if (include.nobs == TRUE) {
    n <- length(model$linear.predictors)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\\ obs.")
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
    include.adjrs = TRUE, include.nobs = TRUE, beside = FALSE,
    include.suffix = FALSE, ...) {
  if (beside == TRUE) {
    equationList <- list()
    for(eq in model$eq) {  # go through estimated equations
      s <- summary(eq, ...)  # extract model summary
      names <- rownames(coef(s))
      co <- coef(s)[, 1]
      se <- coef(s)[, 2]
      pval <- coef(s)[, 4]

      gof <- numeric()
      gof.names <- character()
      gof.decimal <- logical()
      if (include.rsquared == TRUE) {
        rs <- s$r.squared  # extract r-squared
        gof <- c(gof, rs)
        gof.names <- c(gof.names, "R$^2$")
        gof.decimal <- c(gof.decimal, TRUE)
      }
      if (include.adjrs == TRUE) {
        adj <- s$adj.r.squared  # extract adjusted r-squared
        gof <- c(gof, adj)
        gof.names <- c(gof.names, "Adj.\ R$^2$")
        gof.decimal <- c(gof.decimal, TRUE)
      }
      if (include.nobs == TRUE) {
        n <- length(s$residuals)  # extract number of observations
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
        gof.decimal = gof.decimal,
        model.name = s$eqnLabel
      )
      equationList[[eq$eqnNo]] <- tr
    }
    return(equationList)  #returns a list of table.content lists
  } else {
    nm <- character()
    co <- se <- pval <- numeric()
    for(eq in model$eq) {  # go through estimated equations
      s <- summary(eq, ...)  # extract model summary
      nm <- c(nm, sapply(rownames(coef(s)), function(x) {
        if (include.suffix == FALSE) {
          paste0(eq$eqnLabel, ": ", x)
        } else {
          paste0(x, " (", eq$eqnLabel, ")")
        }
      }))
      co <- c(co, coef(s)[, 1])
      se <- c(se, coef(s)[, 2])
      pval <- c(pval, coef(s)[, 4])
    }

    gof <- numeric()
    gof.names <- character()
    gof.decimal <- logical()
    for(eq in model$eq) {
      s <- summary(eq, ...)
      if (include.rsquared == TRUE) {
        rs <- s$r.squared  # extract r-squared
        gof <- c(gof, rs)
        gof.names <- c(gof.names, ifelse(include.suffix == TRUE,
            paste0("R$^2$ (", eq$eqnLabel, ")"),
            paste0(eq$eqnLabel, ": R$^2$")))
        gof.decimal <- c(gof.decimal, TRUE)
      }
    }

    for(eq in model$eq) {
      s <- summary(eq, ...)
      if (include.adjrs == TRUE) {
        adj <- s$adj.r.squared  # extract adjusted r-squared
        gof <- c(gof, adj)
        gof.names <- c(gof.names, ifelse(include.suffix == TRUE,
            paste0("Adj.\ R$^2$ (", eq$eqnLabel, ")"),
            paste0(eq$eqnLabel, ": Adj.\ R$^2$")))
        gof.decimal <- c(gof.decimal, TRUE)
      }
    }

    if (include.nobs == TRUE) {  # number of observations
      gof <- c(gof, nobs(model))
      gof.names <- c(gof.names, "Num.\ obs.\ (total)")
      gof.decimal <- c(gof.decimal, FALSE)
    }

    tr <- createTexreg(
      coef.names = nm,
      coef = co,
      se = se,
      pvalues = pval,
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal
    )
    return(tr)
  }
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
    n <- nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.censnobs == TRUE) {
    censnobs <- s$n
    censnobs.names <- names(censnobs)
    gof <- c(gof, censnobs)
    gof.names <- c(gof.names, censnobs.names)
    gof.decimal <- c(gof.decimal, rep(FALSE, length(censnobs)))
  }
  if (include.wald == TRUE) {
    wald <- s$wald
    gof <- c(gof, wald)
    gof.names <- c(gof.names, "Wald Test")
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

setMethod("extract", signature = className("tobit", "AER"),
    definition = extract.tobit)


# extension for vglm objects (VGAM package)
# please report errors to Christoph Riedl at Northeastern University;
# e-mail: c.riedl@neu.edu
extract.vglm <- function(model, include.loglik = TRUE, include.df = TRUE,
    include.nobs = TRUE, ...) {

  s <- VGAM::summary(model)
	names <- rownames(VGAM::coef(s))
	co <- s@coef3[, 1]
	se <- s@coef3[, 2]
	pval <- s@coef3[, 4]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.loglik == TRUE) {
    gof <- c(gof, VGAM::logLik.vlm(model))
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.df == TRUE) {
    gof <- c(gof, df <- s@df[2])
    gof.names <- c(gof.names, "DF")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, length(VGAM::residuals(s)))
    gof.names <- c(gof.names, "Num.\\ obs.")
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

setMethod("extract", signature = className("vglm", "VGAM"),
    definition = extract.vglm)


# extension for weibreg objects (eha package)
extract.weibreg <- function(model, include.aic = TRUE, include.loglik = TRUE,
    include.lr = TRUE, include.nobs = TRUE, include.events = TRUE,
    include.trisk = TRUE, ...) {

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
  if (include.aic == TRUE) {
    aic <- 2 * model$loglik[2] + 2 * length(coefs)
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, FALSE)
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

extract.coxreg <- extract.weibreg
setMethod("extract", signature = className("coxreg", "eha"),
    definition = extract.coxreg)


# extension for wls objects (metaSEM package)
# please report errors to Christoph Riedl at Northeastern University;
# e-mail: c.riedl@neu.edu
extract.wls <- function(model, include.nobs = TRUE, ...) {

	coefnames <- rownames(summary(model)$coef)
	coefs <- summary(model)$coef[, 1]
	se <- as.numeric(summary(model)$coef[, 2])
	pval <- summary(model)$coef[, 6]

  # Compute average variance extracted
	# Based on: http://openmx.psyc.virginia.edu/thread/3988
	# Could also check description of reliability() from {semTools}
	mat <- model$mx.fit$impliedS1$result
	if (is.null(mat)) {
	  ave <- NULL
	} else {
  	ave <- mean(mat[nrow(mat), -ncol(mat)])
	}

	chi      <- summary(model)$stat["Chi-square of independence model", 1]
	dfs       <- summary(model)$stat["DF of independence model", 1]
	# chi.pval <- summary(model)$stat["p value of target model", 1]
	# if(pval < .0001) pval <- "< .0001"
	rmsea    <- summary(model)$stat["RMSEA", 1]
	rmseall  <- summary(model)$stat["RMSEA lower 95% CI", 1]
	rmseaul  <- summary(model)$stat["RMSEA upper 95% CI", 1]
	cfi      <- summary(model)$stat["CFI", 1]

	gof <- c(chi, dfs, rmsea, rmseall, rmseaul, cfi)
	gof.names <- c("Chi-square of independence model",
	    "DF of independence model", "RMSEA", "RMSEA lower 95 percent CI",
	    "RMSEA upper 95 percent CI", "CFI")
  gof.decimal <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
  if (!is.null(ave)) {
    gof <- c(gof, ave)
    gof.names <- c(gof.names, "Average variance extracted")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, summary(model)$stat["Sample size", 1])
    gof.names <- c(gof.names, "Num.\\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }

	tr <- createTexreg(
		coef.names = coefnames,
		coef = coefs,
		se = se,
		pvalues = pval,
		gof.names = gof.names,
		gof = gof,
		gof.decimal = gof.decimal
	)
	return(tr)
}

setMethod("extract", signature = className("wls", "metaSEM"),
    definition = extract.wls)


# extension for zelig objects (Zelig package < 5.0)
extract.zelig <- function(model, include.aic = TRUE, include.bic = TRUE,
    include.loglik = TRUE, include.deviance = TRUE, include.nobs = TRUE,
    include.rsquared = TRUE, include.adjrs = TRUE, include.fstatistic = TRUE,
    ...) {

  s <- summary(model, ...)

  if ("relogit" %in% class(model) || "logit" %in% class(model) ||
      "ls" %in% class(model) || "probit" %in% class(model) ||
      "ologit" %in% class(model)) {
    coefficient.names <- rownames(s$coef)
    coefficients <- s$coef[, 1]
    standard.errors <- s$coef[, 2]
    if ("ologit" %in% class(model)) {
      tval <- s$coef[, 3]
      significance <- 2 * pt(-abs(tval), s$df.residual)
    } else {
      significance <- s$coef[, 4]
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
      gof.names <- c(gof.names, "Log Likelihood")
      gof.decimal <- c(gof.decimal, TRUE)
    }
    if (include.deviance == TRUE) {
      dev <- s$deviance
      if (!is.null(dev)) {
        gof <- c(gof, dev)
        gof.names <- c(gof.names, "Deviance")
        gof.decimal <- c(gof.decimal, TRUE)
      }
    }
    if (include.nobs == TRUE) {
      n <- nrow(model$data)
      gof <- c(gof, n)
      gof.names <- c(gof.names, "Num.\ obs.")
      gof.decimal <- c(gof.decimal, FALSE)
    }
    if (include.rsquared == TRUE) {
      rs <- s$r.squared  #extract R-squared
      if (!is.null(rs)) {
        gof <- c(gof, rs)
        gof.names <- c(gof.names, "R$^2$")
        gof.decimal <- c(gof.decimal, TRUE)
      }
    }
    if (include.adjrs == TRUE) {
      adj <- s$adj.r.squared  #extract adjusted R-squared
      if (!is.null(adj)) {
        gof <- c(gof, adj)
        gof.names <- c(gof.names, "Adj.\ R$^2$")
        gof.decimal <- c(gof.decimal, TRUE)
      }
    }
    if (include.fstatistic == TRUE) {
      fstat <- s$fstatistic[[1]]
      if (!is.null(fstat)) {
        gof <- c(gof, fstat)
        gof.names <- c(gof.names, "F statistic")
        gof.decimal <- c(gof.decimal, TRUE)
      }
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
  } else if ("mlogit" %in% class(model)) {
    coefficient.names <- rownames(s@coef3)
    coefficients <- s@coef3[, 1]
    standard.errors <- s@coef3[, 2]
    zval <- s@coef3[, 3]
    significance <- 2 * pnorm(abs(zval), lower.tail = FALSE)

    gof <- numeric()
    gof.names <- character()
    gof.decimal <- logical()
    if (include.loglik == TRUE) {
      lik <- logLik(model)[1]
      gof <- c(gof, lik)
      gof.names <- c(gof.names, "Log Likelihood")
      gof.decimal <- c(gof.decimal, TRUE)
    }
    if (include.deviance == TRUE) {
      dev <- deviance(s)
      if (!is.null(dev)) {
        gof <- c(gof, dev)
        gof.names <- c(gof.names, "Deviance")
        gof.decimal <- c(gof.decimal, TRUE)
      }
    }
    if (include.nobs == TRUE) {
      n <- nrow(model$data)
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
  } else if ("tobit" %in% class(model)) {
    coefficient.names <- rownames(s$table)
    coefficients <- s$table[, 1]
    standard.errors <- s$table[, 2]
    significance <- s$table[, 5]
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
    if (include.nobs == TRUE) {
        n <- nrow(model$data)
        gof <- c(gof, n)
        gof.names <- c(gof.names, "Num.\ obs.")
        gof.decimal <- c(gof.decimal, FALSE)
    }
    tr <- createTexreg(coef.names = coefficient.names, coef = coefficients,
        se = standard.errors, pvalues = significance, gof.names = gof.names,
        gof = gof, gof.decimal = gof.decimal)
    return(tr)
  } else {
    stop(paste("Only the following Zelig models are currently supported:",
        "logit, ls, mlogit, ologit, probit, relogit, tobit."))
  }
}

setMethod("extract", signature = className("zelig", "Zelig"),
    definition = extract.zelig)


# extension for Zelig objects (Zelig package >= 5.0)
extract.Zelig <- function(model, include.nobs = TRUE, include.nimp = TRUE, ...) {
  if (model$mi) {
    if (!exists("combine_coef_se", where = "package:Zelig",
        mode = "function")) {
      stop("texreg relies on Zelig's combine_coef_se function to extract model information. Install Zelig >= 5.0-17 to see if texreg can format your model.")
    }
    combined <- Zelig::combine_coef_se(model, messages = FALSE)
    gof <- gof.names <- gof.decimal <- NULL
    if (include.nobs) {
      gof <- c(gof, nrow(model$data))
      gof.names <- c(gof.names, 'Num. obs.')
      gof.decimal <- c(gof.decimal, FALSE)
    }
    if (include.nimp) {
      if (class(model$originaldata)[1] == 'amelia') {
        gof <- c(gof, model$originaldata$m)
        gof.names <- c(gof.names, 'Num. imp.')
        gof.decimal <- c(gof.decimal, FALSE)
      } else if (class(model$originaldata)[1] == 'mi') { # when imputed dataset was created using to_zelig_mi
        gof <- c(gof, length(model$originaldata))
        gof.names <- c(gof.names, 'Num. imp.')
        gof.decimal <- c(gof.decimal, FALSE)
      }
    }
    out <- createTexreg(coef.names = row.names(combined),
                        coef = combined[, 'Estimate'],
                        se = combined[, 'Std.Error'],
                        pvalues = combined[, 'Pr(>|z|)'],
                        gof.names = gof.names,
                        gof = gof,
                        gof.decimal = gof.decimal)
  } else {
    if ("Zelig-relogit" %in% class(model)) { # remove when users update to Zelig 5.0-16
      mod_original <- model$zelig.out$z.out[[1]]
      class(mod_original) <- "glm"
    }
    else if ("Zelig-tobit" %in% class(model)) { # remove when users update to Zelig 5.0-16
      mod_original <- model$zelig.out$z.out[[1]]
    } else {
      if (!exists("from_zelig_model", where = "package:Zelig",
          mode = "function")) {
        stop("texreg relies on Zelig's from_zelig_model function to extract model information. Install Zelig >= 5.0-16 to see if texreg can format your model.")
      }
      mod_original <- try(Zelig::from_zelig_model(model), silent = TRUE)
      if (class(mod_original)[1] == "try-error") {
        stop("texreg relies on Zelig's from_zelig_model function to extract information from Zelig models. from_zelig_model does not appear to support models of class ",
               class(model)[1], ".")
      }
    }
    out <- extract(mod_original, include.nobs = include.nobs, ...)
  }
  return(out)
}

setMethod("extract", signature = className("Zelig", "Zelig"),
    definition = extract.Zelig)


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
          gof.decimal = gof.decimal,
          model.name = "Count model"
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
          gof.decimal = gof.decimal,
          model.name = "Zero model"
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

