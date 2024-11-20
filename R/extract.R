# Generic function -------------------------------------------------------------

#' Extract details from statistical models for table construction
#'
#' Extract details from statistical models for table construction. The function
#' has methods for a range of statistical models.
#'
#' The \code{\link{extract}} function serves to retrieve coefficients, standard
#' errors, p-values, confidence intervals, and goodness-of-fit statistics from
#' statistical models in \R. More than 100 \code{\link{extract}} methods
#' ("extensions") for various statistical models are available. The function
#' returns a \linkS4class{texreg} object.
#'
#' \code{\link{extract}} is a generic function, which extracts coefficients and
#' GOF measures from statistical model objects. \code{\link{extract}} methods
#' for the specific model types are called by the generic \code{\link{extract}}
#' function if it encounters a model known to be handled by the specific method.
#' The output is a \linkS4class{texreg} object, which is subsequently used by
#' the \code{\link{texreg}} function and related functions.
#'
#' To list the model classes for which extract methods exist, type
#' \code{showMethods("extract")} or \code{methods("extract")}. To show the
#' method definition (i.e., the function body) of a specific extract method, use
#' the \code{getMethod} function, for example \code{getMethod("extract", "lm")}
#' for linear models. To get help on a specific extract method, type something
#' like \code{?`extract,lm-method`} (or something similar for other models,
#' where \code{"lm"} needs to be replaced by the class name of the respective
#' model). You can also list the available methods by displaying the
#' \link[=texreg-package]{texreg package} help index.
#'
#' Users can contribute their own extensions for additional statistical models.
#' Examples are contained in the article in the Journal of Statistical Software
#' referenced below. Suggestions can be submitted as pull requests at
#' \url{https://github.com/leifeld/texreg/pulls} or through the issue tracker at
#' \url{https://github.com/leifeld/texreg/issues}.
#'
#' @param model A statistical model object.
#' @param ... Custom parameters, which are handed over to subroutines. The
#'   arguments are usually passed to the \code{summary} function, but in some
#'   cases to other functions.
#' @return The function returns a \linkS4class{texreg} object.
#'
#' @author Philip Leifeld
#' @seealso \code{\link{createTexreg}}, \code{\link{matrixreg}},
#'   \code{\link{screenreg}}, \code{\link{texreg}}
#'
#' @references Leifeld, Philip (2013). texreg: Conversion of Statistical Model
#'   Output in R to LaTeX and HTML Tables. Journal of Statistical Software
#'   55(8): 1-24. \doi{10.18637/jss.v055.i08}.
#'
#' @export
setGeneric("extract", function(model, ...) standardGeneric("extract"),
           package = "texreg")


# -- extract.Arima (stats) -----------------------------------------------------

#' @noRd
extract.Arima <- function(model,
                          include.pvalues = TRUE,
                          include.aic = TRUE,
                          include.bic = TRUE,
                          include.loglik = TRUE,
                          include.nobs = TRUE,
                          ...) {

  mask <- model$mask
  nam <- names(model$coef)
  co <- model$coef
  se <- sqrt(diag(model$var.coef))

  if (include.pvalues == TRUE) {
    t.rat <- rep(NA, length(mask))
    t.rat[mask] <- co[mask] / se
    p <- 2 * pnorm(-abs(t.rat))
    setmp <- rep(NA, length(mask))
    setmp[mask] <- se
  } else {
    p <- numeric()
    setmp <- se
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
    gof <- c(gof, BIC(model))
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
    gof <- c(gof, model$nobs)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }

  tr <- createTexreg(
    coef.names = nam,
    coef = co,
    se = setmp,
    pvalues = p,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
  return(tr)
}

#' \code{\link{extract}} method for \code{Arima} objects
#'
#' \code{\link{extract}} method for \code{Arima} objects created by the
#' \code{\link[stats]{arima}} function in the \pkg{stats} package.
#'
#' @param model A statistical model object.
#' @param include.pvalues Report p-values?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract Arima
#' @aliases extract.Arima
#' @importFrom stats pnorm BIC
#' @export
setMethod("extract", signature = className("Arima", "stats"),
          definition = extract.Arima)


# -- extract.forecast_ARIMA (forecast) -----------------------------------------

#' @noRd
extract.forecast_ARIMA <- function (model,
                                    include.pvalues = TRUE,
                                    include.aic = TRUE,
                                    include.aicc = TRUE,
                                    include.bic = TRUE,
                                    include.loglik = TRUE,
                                    include.nobs = TRUE,
                                    ...) {
  mask <- model$mask
  nam <- names(model$coef)
  co <- model$coef
  se <- sqrt(diag(model$var.coef))
  if (include.pvalues == TRUE) {
    t.rat <- rep(NA, length(mask))
    t.rat[mask] <- co[mask] / se
    p <- 2 * pnorm(-abs(t.rat))
    setmp <- rep(NA, length(mask))
    setmp[mask] <- se
  } else {
    p <- numeric()
    setmp <- se
  }

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    aic <- model$aic
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
  if (include.nobs == TRUE) {
    gof <- c(gof, model$nobs)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  tr <- createTexreg(
    coef.names = nam,
    coef = co,
    se = setmp,
    pvalues = p,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
  return(tr)
}

#' \code{\link{extract}} method for \code{forecast_ARIMA} objects
#'
#' \code{\link{extract}} method for \code{forecast_ARIMA} objects created by the
#' \code{\link[forecast]{Arima}} function in the \pkg{forecast} package.
#'
#' @param model A statistical model object.
#' @param include.pvalues Report p-values?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.aicc Report AICC in the GOF block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract forecast_ARIMA
#' @aliases extract.forecast_ARIMA
#' @importFrom stats pnorm
#' @export
setMethod("extract",
          signature = className("forecast_ARIMA", "forecast"),
          definition = extract.forecast_ARIMA)


# -- extract.averaging (MuMIn) -------------------------------------------------

#' @noRd
extract.averaging <- function(model, use.ci = FALSE, adjusted.se = FALSE,
                              include.nobs = TRUE, ...) {

  # MuMIn >= 1.15.0 : c("coefmat.subset", "coefmat.full")
  # MuMIn < 1.15.0 : c("coefmat", "coefmat.full")
  ct <- summary(model)$coefmat.full
  coefs <- ct[, 1L]
  se <- ct[, if(adjusted.se) 3L else 2L]

  if (include.nobs == TRUE) {
    gof <- as.numeric(attr(model, "nobs"))
    gof.names <- "Num. obs."
    gof.decimal <- FALSE
  } else {
    gof <- numeric(0L)
    gof.names <- character(0L)
    gof.decimal <- logical(0L)
  }

  if (use.ci == TRUE) {
    ci <- stats::confint(model, full = TRUE)
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

#' \code{\link{extract}} method for \code{averaging} objects
#'
#' \code{\link{extract}} method for \code{averaging} objects created by the
#' \code{\link[MuMIn]{model.avg}} function in the \pkg{MuMIn} package.
#'
#' @param model A statistical model object.
#' @param use.ci Report confidence intervals in the GOF block?
#' @param adjusted.se Report adjusted standard error in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract averaging
#' @aliases extract.averaging
#' @importFrom stats confint
#' @export
setMethod("extract", signature = className("averaging", "MuMIn"),
          definition = extract.averaging)


# -- extract.bergm (Bergm) -----------------------------------------------------

#' @noRd
extract.bergm <- function(model, posterior.median = FALSE, level = 0.95, ...) {
  rNames <- paste("theta",
                  seq(1, model$dim),
                  " (",
                  model$specs[seq(1, model$dim)],
                  ")",
                  sep = "") # from Bergm 5.0.5 summary.Bergm method
  coefs <- apply(model$Theta, 2, ifelse(isTRUE(posterior.median), median, mean))
  alpha <- (1 - level) / 2
  cil <- apply(model$Theta, 2, function(x) quantile(x, (1 - level) / 2))
  ciu <- apply(model$Theta, 2, function(x) quantile(x, 1 - ((1 - level) / 2)))

  tr <- createTexreg(
    coef.names = rNames,
    coef = coefs,
    ci.low = cil,
    ci.up = ciu
  )
  return(tr)
}

#' \code{\link{extract}} method for \code{bergm} objects
#'
#' \code{\link{extract}} method for \code{bergm} objects created by the
#' \code{\link[Bergm]{bergm}} function in the \pkg{Bergm} package.
#'
#' @param model A statistical model object.
#' @param posterior.median Report the posterior median instead of the default
#'   posterior mean as coefficients?
#' @param level Confidence level, i.e., the proportion of the posterior
#'   distribution to be included in the credible interval.
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract bergm
#' @aliases extract.bergm
#' @export
setMethod("extract", signature = className("bergm", "Bergm"),
          definition = extract.bergm)


# -- extract.betamfx (mfx) -----------------------------------------------------

#' @noRd
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{betamfx} objects
#'
#' \code{\link{extract}} method for \code{betamfx} objects created by the
#' \code{\link[mfx]{betamfx}} function in the \pkg{mfx} package.
#'
#' @param model A statistical model object.
#' @param include.pseudors Report pseudo R^2 in the GOF block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract betamfx
#' @aliases extract.betamfx
#' @export
setMethod("extract", signature = className("betamfx", "mfx"),
          definition = extract.betamfx)


# -- extract.betaor (mfx) -----------------------------------------------------

#' @noRd
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{betaor} objects
#'
#' \code{\link{extract}} method for \code{betaor} objects created by the
#' \code{\link[mfx]{betaor}} function in the \pkg{mfx} package.
#'
#' @param model A statistical model object.
#' @param include.pseudors Report pseudo R^2 in the GOF block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract betaor
#' @aliases extract.betaor
#' @export
setMethod("extract", signature = className("betaor", "mfx"),
          definition = extract.betaor)


# -- extract.bife (bife) ----------------------------------------------------

#' @noRd
extract.bife <- function(model,
                         include.loglik = TRUE,
                         include.deviance = TRUE,
                         include.nobs = TRUE,
                         ...) {
  s <- summary(model)
  coefficient.names <- rownames(s$cm)
  co <- s$cm[, 1]
  se <- s$cm[, 2]
  pval <- s$cm[, 4]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.loglik == TRUE) {
    lik	<- logLik(model)
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    gof <- c(gof, deviance(model))
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    n <- s$nobs["nobs"]
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }

  tr <- createTexreg(
    coef.names = coefficient.names,
    coef = co,
    se = se,
    pvalues = pval,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
  return(tr)
}

#' \code{\link{extract}} method for \code{bife} objects
#'
#' \code{\link{extract}} method for \code{bife} objects created by the
#' \code{\link[bife]{bife}} function in the \pkg{bife} package.
#'
#' @param model A statistical model object.
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.deviance Report the residual deviance?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract bife
#' @aliases extract.bife
#' @author Philip Leifeld, Christoph Riedl, Claudia Zucca
#' @importFrom stats logLik
#' @export
setMethod("extract", signature = className("bife", "bife"),
          definition = extract.bife)


# -- extract.broom (broom) -----------------------------------------------------

#' @noRd
extract.broom <- function(model, ...) {
  if (!requireNamespace("broom", quietly = TRUE)) {
    stop("texreg does not directly support models of class ",
         class(model)[1],
         ", but it can sometimes use the ``broom`` package to extract model ",
         "information. Call texreg again after installing the ``broom`` ",
         "package to see if this is possible.")
  }
  coefficients <- try(broom::tidy(model)[, c("term",
                                             "estimate",
                                             "std.error",
                                             "p.value")],
                      silent = TRUE)
  gof <- try({
    # extract
    out <- broom::glance(model)[1, ]
    gof.decimal <- sapply(out, function(k) class(k)[1]) # type inference
    gof.decimal <- ifelse(gof.decimal %in% c("integer", "logical"), FALSE, TRUE)
    out <- data.frame("gof.names" = colnames(out),
                      "gof" = as.numeric(out),
                      "gof.decimal" = gof.decimal,
                      stringsAsFactors = FALSE)
    # rename
    gof_dict <- c(
      "adj.r.squared" = "Adj. R$^2$",
      "deviance" = "Deviance",
      "df" = "DF",
      "df.residual" = "DF Resid.",
      "finTol" = "Tolerance",
      "isConv" = "Convergence",
      "logLik" = "Log Likelihood",
      "null.deviance" = "Deviance (Null)",
      "p.value" = "P Value",
      "r.squared" = "R$^2$",
      "sigma" = "Sigma",
      "statistic" = "Statistic"
    )
    gof_dict <- gof_dict[names(gof_dict) %in% out$gof.names]
    idx <- match(names(gof_dict), out$gof.names)
    out$gof.names[idx] <- gof_dict
    if (any(is.na(out$gof))) {
      warning(paste("texreg used the broom package to extract the following GOF",
                    "measures, but could not cast them to numeric type:",
                    out$gof.names[is.na(out$gof)]))
    }
    stats::na.omit(out)
  }, silent = TRUE)
  if ("try-error" %in% class(coefficients) || "try-error" %in% class(gof)) {
    stop("Neither texreg nor broom supports models of class ",
         class(model)[1],
         ".")
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

#' \code{\link{extract}} method for \code{broom} objects
#'
#' \code{\link{extract}} method for \code{broom} objects created by the
#' \code{\link[broom]{broom}} function in the \pkg{broom} package.
#'
#' @param model A statistical model object.
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract broom
#' @aliases extract.broom extract.ANY extract.ANY-method extract,ANY-method
#' @importFrom stats na.omit
#' @export
setMethod("extract", signature = className("ANY"), definition = extract.broom)


# -- extract.biglm (biglm) -----------------------------------------------------

#' @noRd
extract.biglm <- function(model, include.nobs = TRUE, include.aic = TRUE,
                          use.ci = FALSE, ...) {

  tab <-summary(model)$mat

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    gof <- c(gof, model$n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.aic == TRUE) {
    gof <- c(gof, AIC(model))
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  if (use.ci == TRUE) {
    tr <- createTexreg(
      coef.names = rownames(tab),
      coef = tab[, 1],
      ci.low = tab[, 2],
      ci.up = tab[, 3],
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal
    )

  } else {
    tr <- createTexreg(
      coef.names = rownames(tab),
      coef = tab[, 1],
      se = tab[, 4],
      pvalues = tab[, 5],
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal
    )
  }
  return(tr)
}

#' \code{\link{extract}} method for \code{biglm} objects
#'
#' \code{\link{extract}} method for \code{biglm} objects created by the
#' \code{\link[biglm]{biglm}} function in the \pkg{biglm} package.
#'
#' @param model A statistical model object.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param use.ci Report confidence intervals in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract biglm
#' @aliases extract.biglm
#' @author Claudia Zucca, Philip Leifeld
#' @export
setMethod("extract", signature = className("biglm", "biglm"),
          definition = extract.biglm)


# -- extract.brmsfit (brms) ----------------------------------------------------

#' @noRd
extract.brmsfit <- function (model,
                             use.HDI = TRUE,
                             level = 0.9,
                             include.random = TRUE,
                             include.rsquared = TRUE,
                             include.nobs = TRUE,
                             include.loo.ic = TRUE,
                             reloo = FALSE,
                             include.waic = TRUE,
                             ...) {
  sf <- summary(model, ...)$fixed
  coefnames <- rownames(sf)
  coefs <- sf[, 1]
  se <- sf[, 2]
  if (isTRUE(use.HDI)) {
    hdis <- coda::HPDinterval(brms::as.mcmc(model, combine_chains = TRUE),
                              prob = level)
    hdis <- hdis[seq(1:length(coefnames)), ]
    ci.low = hdis[, "lower"]
    ci.up = hdis[, "upper"]
  } else { # default using 95% posterior quantiles from summary.brmsfit
    ci.low = sf[, 3]
    ci.up = sf[, 4]
  }

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (isTRUE(include.random) & isFALSE(!nrow(model$ranef))) {
    sr <- summary(model, ...)$random
    sd.names <- character()
    sd.values <- numeric()
    for (i in 1:length(sr)) {
      sd <- sr[[i]][, 1]
      sd.names <- c(sd.names, paste0("SD: ", names(sr)[[i]], names(sd)))
      sd.values <- c(sd.values, sd)
    }
    gof <- c(gof, sd.values)
    gof.names <- c(gof.names, sd.names)
    gof.decimal <- c(gof.decimal, rep(TRUE, length(sd.values)))
  }
  if (isTRUE(include.rsquared)) {
    rs <- brms::bayes_R2(model)[1]
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (isTRUE(include.nobs)) {
    n <- stats::nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (isTRUE(include.loo.ic)) {
    looic <- brms::loo(model, reloo = reloo)$estimates["looic", "Estimate"]
    gof <- c(gof, looic)
    gof.names <- c(gof.names, "loo IC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (isTRUE(include.waic)) {
    waic <- brms::waic(model)$estimates["waic", "Estimate"]
    gof <- c(gof, waic)
    gof.names <- c(gof.names, "WAIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  tr <- createTexreg(coef.names = coefnames,
                     coef = coefs,
                     se = se,
                     ci.low = ci.low,
                     ci.up = ci.up,
                     gof.names = gof.names,
                     gof = gof,
                     gof.decimal = gof.decimal)
  return(tr)
}

#' \code{\link{extract}} method for \code{brmsfit} objects
#'
#' \code{\link{extract}} method for \code{brmsfit} objects created by the
#' \code{\link[brms]{brm}} function in the \pkg{brms} package.
#'
#' @param model A statistical model object.
#' @param use.HDI Report highest posterior density (HPD) intervals (HDI) using
#'   the \code{\link[coda]{HPDinterval}} function in the \pkg{coda} package,
#'   with the probability given in the \code{level} argument, instead of the
#'   default 95 percent posterior quantiles?
#' @param level Significance level (\code{1 - alpha}) for HPDs (in combination
#'   with the \code{use.HDI} argument).
#' @param include.random Include random effects (standard deviations) in the GOF
#'   block of the table?
#' @param include.rsquared Report R^2 in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.loo.ic Report Leave-One-Out Information Criterion?
#' @param reloo Recompute exact cross-validation for problematic observations
#'   for which approximate leave-one-out cross-validation may return incorrect
#'   results? This is done using the \code{\link[brms]{reloo.brmsfit}} function
#'   and may take some time to compute.
#' @param include.waic Report Widely Applicable Information Criterion (WAIC)?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract brmsfit
#' @aliases extract.brmsfit
#' @author Hyunjin (Jin) Song, Philip Leifeld
#' @importFrom stats nobs
#' @export
setMethod("extract",
          signature = className("brmsfit", "brms"),
          definition = extract.brmsfit)


# -- extract.btergm (btergm) ---------------------------------------------------

#' @noRd
extract.btergm <- function(model, level = 0.95, include.nobs = TRUE, ...) {

  tab <- btergm::confint(model, level = level)

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    gof <- c(gof, model@nobs)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }

  tr <- createTexreg(
    coef.names = rownames(tab),
    coef = tab[, which(grepl("Estimate", colnames(tab)))],
    ci.low = tab[, which(grepl("2.5", colnames(tab)))],
    ci.up = tab[, which(grepl("97.5", colnames(tab)))],
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )

  return(tr)
}

#' \code{\link{extract}} method for \code{btergm} objects
#'
#' \code{\link{extract}} method for \code{btergm} objects created by the
#' \code{\link[btergm]{btergm}} function in the \pkg{btergm} package.
#'
#' @param model A statistical model object.
#' @param level Significance or confidence level (\code{1 - alpha}) for
#'   computing confidence intervals.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract btergm
#' @aliases extract.btergm
#' @export
setMethod("extract", signature = className("btergm", "btergm"),
          definition = extract.btergm)


# -- extract.censReg (censReg) -------------------------------------------------

#' @noRd
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
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, s$nObs)
    gof.names <- c(gof.names, "Num. obs.", "Left-censored", "Uncensored",
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

#' \code{\link{extract}} method for \code{censReg} objects
#'
#' \code{\link{extract}} method for \code{censReg} objects created by the
#' \code{censReg} function in the \pkg{censReg} package.
#'
#' @param model A statistical model object.
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract censReg
#' @aliases extract.censReg
#' @export
setMethod("extract", signature = className("censReg", "censReg"),
          definition = extract.censReg)


# -- extract.clm (ordinal) -----------------------------------------------------

#' @noRd
extract.clm <- function(model,
                        include.thresholds = TRUE,
                        include.aic = TRUE,
                        include.bic = TRUE,
                        include.loglik = TRUE,
                        include.nobs = TRUE,
                        ...) {
  s <- summary(model, ...)

  tab <- s$coefficients
  thresh <- tab[rownames(tab) %in% names(s$aliased$alpha), ]
  threshold.names <- rownames(thresh)
  threshold.coef <- thresh[, 1]
  threshold.se <- thresh[, 2]
  threshold.pval <- thresh[, 4]
  beta <- tab[rownames(tab) %in% names(s$aliased$beta), , drop = FALSE]
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{clm} objects
#'
#' \code{\link{extract}} method for \code{clm} objects created by the
#' \code{\link[ordinal]{clm}} function in the \pkg{ordinal} package.
#'
#' @param model A statistical model object.
#' @param include.thresholds Report thresholds in the GOF block?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract clm
#' @aliases extract.clm
#' @export
setMethod("extract", signature = className("clm", "ordinal"),
          definition = extract.clm)


# -- extract.clmm (ordinal) ----------------------------------------------------

#' @noRd
extract.clmm <- function(model,
                         include.thresholds = TRUE,
                         include.loglik = TRUE,
                         include.aic = TRUE,
                         include.bic = TRUE,
                         include.nobs = TRUE,
                         include.groups = TRUE,
                         include.variance = TRUE,
                         ...) {
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{clmm} objects
#'
#' \code{\link{extract}} method for \code{clmm} objects created by the
#' \code{\link[ordinal]{clmm}} function in the \pkg{ordinal} package.
#'
#' @param model A statistical model object.
#' @param include.thresholds Report thresholds in the GOF block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.groups Report the number of groups?
#' @param include.variance Report group variances?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract clmm
#' @aliases extract.clmm
#' @export
setMethod("extract", signature = className("clmm", "ordinal"),
          definition = extract.clmm)


# -- extract.clogit (survival) -------------------------------------------------

#' @noRd
extract.clogit <- function(model,
                           include.aic = TRUE,
                           include.rsquared = TRUE,
                           include.maxrs = TRUE,
                           include.events = TRUE,
                           include.nobs = TRUE,
                           include.missings = TRUE,
                           ...) {
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

  event <- model$nevent
  n <- model$n
  mis <- length(model$na.action)
  rs <- s$rsq[1]
  maxrs <- s$rsq[2]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    aic <- stats::extractAIC(model)[2]
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
    gof.names <- c(gof.names, "Max. R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.events == TRUE) {
    gof <- c(gof, event)
    gof.names <- c(gof.names, "Num. events")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{clogit} objects
#'
#' \code{\link{extract}} method for \code{clogit} objects created by the
#' \code{\link[survival]{clogit}} function in the \pkg{survival} package.
#'
#' @param model A statistical model object.
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.rsquared Report R^2 in the GOF block?
#' @param include.maxrs Report maximal R^2 in the GOF block?
#' @param include.events Report the number of events in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.missings Report number of missing data points in the GOF
#'   block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract clogit
#' @aliases extract.clogit
#' @importFrom stats extractAIC
#' @export
setMethod("extract", signature = className("clogit", "survival"),
          definition = extract.clogit)


# -- extract.coeftest (lmtest) -------------------------------------------------

#' @noRd
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

#' \code{\link{extract}} method for \code{coeftest} objects
#'
#' \code{\link{extract}} method for \code{coeftest} objects created by the
#' \code{\link[lmtest]{coeftest}} function in the \pkg{lmtest} package.
#'
#' @param model A statistical model object.
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract coeftest
#' @aliases extract.coeftest
#' @export
setMethod("extract", signature = className("coeftest", "lmtest"),
          definition = extract.coeftest)


# -- extract.coxph (survival) -------------------------------------------------

#' @noRd
extract.coxph <- function(model,
                          include.aic = TRUE,
                          include.rsquared = TRUE,
                          include.maxrs = TRUE,
                          include.events = TRUE,
                          include.nobs = TRUE,
                          include.missings = TRUE,
                          include.zph = TRUE,
                          ...) {
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

  event <- model$nevent
  n <- model$n
  mis <- length(model$na.action)
  rs <- s$rsq[1]
  maxrs <- s$rsq[2]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    aic <- stats::extractAIC(model)[2]
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
    gof.names <- c(gof.names, "Max. R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.events == TRUE) {
    gof <- c(gof, event)
    gof.names <- c(gof.names, "Num. events")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{coxph} objects
#'
#' \code{\link{extract}} method for \code{coxph} objects created by the
#' \code{\link[survival]{coxph}} function in the \pkg{survival} package.
#'
#' @param model A statistical model object.
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.rsquared Report R^2 in the GOF block?
#' @param include.maxrs Report maximal R^2 in the GOF block?
#' @param include.events Report the number of events in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.missings Report number of missing data points in the GOF
#'   block?
#' @param include.zph Report proportional hazard test in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract coxph
#' @aliases extract.coxph
#' @importFrom stats extractAIC pchisq
#' @export
setMethod("extract", signature = className("coxph", "survival"),
          definition = extract.coxph)


# -- extract.coxph.penal (survival) --------------------------------------------

#' @noRd
extract.coxph.penal <- function(model,
                                include.aic = TRUE,
                                include.rsquared = TRUE,
                                include.maxrs = TRUE,
                                include.events = TRUE,
                                include.nobs = TRUE,
                                include.missings = TRUE,
                                include.zph = TRUE,
                                ...) {

  coefficients <- coef(model, ...)
  coefficient.names <- names(coefficients)
  standard.errors <- sqrt(diag(model$var))
  significance <- 1 - pchisq(model$coefficients^2 / diag(model$var), 1)

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
    aic <- stats::extractAIC(model)[2]
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
    gof.names <- c(gof.names, "Max. R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.events == TRUE) {
    gof <- c(gof, event)
    gof.names <- c(gof.names, "Num. events")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{coxph.penal} objects
#'
#' \code{\link{extract}} method for \code{coxph.penal} objects created by the
#' \code{\link[survival]{coxph}} function in the \pkg{survival} package.
#'
#' @inheritParams extract,coxph-method
#'
#' @method extract coxph.penal
#' @aliases extract.coxph.penal
#' @export
setMethod("extract", signature = className("coxph.penal", "survival"),
          definition = extract.coxph.penal)


# -- extract.ergm (ergm) -------------------------------------------------------

#' @noRd
extract.ergm <- function(model, include.aic = TRUE, include.bic = TRUE,
                         include.loglik = TRUE, ...) {
  s <- summary(model, ...)
  coefs <- coef(s)
  coefficient.names <- rownames(coefs)
  coefficients <- coefs[, 1]
  standard.errors <- coefs[, 2]
  significance <- coefs[, 5]

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

#' \code{\link{extract}} method for \code{ergm} objects
#'
#' \code{\link{extract}} method for \code{ergm} objects created by the
#' \code{\link[ergm]{ergm}} function in the \pkg{ergm} package.
#'
#' @param model A statistical model object.
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract ergm
#' @aliases extract.ergm
#' @export
setMethod("extract", signature = className("ergm", "ergm"),
          definition = extract.ergm)


# -- extract.ergmm (latentnet) -------------------------------------------------

#' @noRd
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

#' \code{\link{extract}} method for \code{ergmm} objects
#'
#' \code{\link{extract}} method for \code{ergmm} objects created by the
#' \code{\link[latentnet]{ergmm}} function in the \pkg{latentnet}
#' package.
#'
#' @param model A statistical model object.
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract ergmm
#' @aliases extract.ergmm
#' @export
setMethod("extract", signature = className("ergmm", "latentnet"),
          definition = extract.ergmm)


# -- extract.ets (forecast) ----------------------------------------------------

#' @noRd
extract.ets <- function (model,
                         include.pvalues = FALSE,
                         include.aic = TRUE,
                         include.aicc = TRUE,
                         include.bic = TRUE,
                         include.loglik = TRUE,
                         ...) {
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

#' \code{\link{extract}} method for \code{ets} objects
#'
#' \code{\link{extract}} method for \code{ets} objects created by the
#' \code{\link[forecast]{ets}} function in the \pkg{forecast} package.
#'
#' @param model A statistical model object.
#' @param include.pvalues Report p-values?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.aicc Report AICC in the GOF block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract ets
#' @aliases extract.ets
#' @export
setMethod("extract",
          signature = className("ets", "forecast"),
          definition = extract.ets)


# -- extract.fGARCH (fGarch) ---------------------------------------------------

#' @noRd
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{fGARCH} objects
#'
#' \code{\link{extract}} method for \code{fGARCH} objects created by the
#' \code{\link[fGarch]{garchFit}} function in the \pkg{fGarch} package.
#'
#' @param model A statistical model object.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract fGARCH
#' @aliases extract.fGARCH
#' @export
setMethod("extract", signature = className("fGARCH", "fGarch"),
          definition = extract.fGARCH)


# -- extract.feis (feisr) ------------------------------------------------------

#' @noRd
extract.feis <- function(model,
                         include.rsquared = TRUE,
                         include.adjrs = TRUE,
                         include.nobs = TRUE,
                         include.groups = TRUE,
                         include.rmse = TRUE,
                         ...) {
  s <- summary(model, ...)

  coefficient.names <- rownames(coef(s))
  coefficients <- coef(s)[, 1]
  standard.errors <- coef(s)[, 2]
  significance <- coef(s)[, 4]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (isTRUE(include.rsquared)) {
    rs <- s$r.squared[1]
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (isTRUE(include.adjrs)) {
    adj <- s$r.squared[2]
    gof <- c(gof, adj)
    gof.names <- c(gof.names, "Adj. R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (isTRUE(include.nobs)) {
    n <- length(model$residuals)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (isTRUE(include.groups)) {
    grps <-length(unique(model$id))
    grp.names <- model$call[[match(c("id"), names(model$call))]]
    grp.names <- paste("Num. groups:", grp.names)
    gof <- c(gof, grps)
    gof.names <- c(gof.names, grp.names)
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (isTRUE(include.rmse)) {
    rmse <- sqrt(sum((model$residuals * model$residuals)) / model$df.residual)
    gof <- c(gof, rmse)
    gof.names <- c(gof.names, "RMSE")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  tr <- texreg::createTexreg(
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

#' \code{\link{extract}} method for \code{feis} objects
#'
#' \code{\link{extract}} method for \code{feis} objects created by the
#' \code{\link[feisr]{feis}} function in the \pkg{feisr} package.
#'
#' @param model A statistical model object.
#' @param include.rsquared Report R^2 in the GOF block?
#' @param include.adjrs Report adjusted R^2 in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.groups Report the number of groups?
#' @param include.rmse Report the root mean square error (RMSE; = residual
#'   standard deviation) in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract feis
#' @aliases extract.feis
#' @author Tobias Rüttenauer, Philip Leifeld
#' @importFrom stats nobs coef
#' @export
setMethod("extract", signature = className("feis", "feisr"),
          definition = extract.feis)


# -- extract.betareg (betareg) -------------------------------------------------

#' @noRd
extract.betareg <- function(model, include.precision = TRUE,
                            include.pseudors = TRUE, include.loglik = TRUE,
                            include.nobs = TRUE, ...) {

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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{betareg} objects
#'
#' \code{\link{extract}} method for \code{betareg} objects created by the
#' \code{\link[betareg]{betareg}} function in the \pkg{betareg} package.
#'
#' @param model A statistical model object.
#' @param include.precision Report precision in the GOF block?
#' @param include.pseudors Report pseudo R^2 in the GOF block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract betareg
#' @aliases extract.betareg
#' @export
setMethod("extract", signature = className("betareg", "betareg"),
          definition = extract.betareg)


# -- extract.felm (lfe) -----------------------------------------------

#' @noRd
extract.felm <- function(model,
                         include.nobs = TRUE,
                         include.rsquared = TRUE,
                         include.adjrs = TRUE,
                         include.fstatistic = FALSE,
                         include.proj.stats = TRUE,
                         include.groups = TRUE,
                         ...) {

  s <- summary(model, ...)
  nam <- rownames(s$coefficients)
  co <- s$coefficients[, 1]
  se <- s$coefficients[, 2]
  pval <- s$coefficients[, 4]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    gof <- c(gof, s$N)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.rsquared == TRUE) {
    gof <- c(gof, s$r2)
    gof.decimal <- c(gof.decimal, TRUE)
    if (include.proj.stats == TRUE) {
      gof <- c(gof, s$P.r.squared)
      gof.decimal <- c(gof.decimal, TRUE)
      gof.names <- c(gof.names, "R$^2$ (full model)", "R$^2$ (proj model)")
    } else {
      gof.names <- c(gof.names, "R$^2$")
    }
  }
  if (include.adjrs == TRUE) {
    gof <- c(gof, s$r2adj)
    gof.decimal <- c(gof.decimal, TRUE)
    if (include.proj.stats == TRUE) {
      gof <- c(gof, s$P.adj.r.squared)
      gof.decimal <- c(gof.decimal, TRUE)
      gof.names <- c(gof.names, "Adj. R$^2$ (full model)", "Adj. R$^2$ (proj model)")
    } else {
      gof.names <- c(gof.names, "Adj. R$^2$")
    }
  }
  if (include.fstatistic == TRUE) {
    gof <- c(gof, s$F.fstat[1], s$F.fstat[4])
    gof.decimal <- c(gof.decimal, TRUE, TRUE)
    if (include.proj.stats == TRUE) {
      gof <- c(gof, s$P.fstat[length(s$P.fstat) - 1], s$P.fstat[1])
      gof.decimal <- c(gof.decimal, TRUE, TRUE)
      gof.names <- c(gof.names, "F statistic (full model)",
                     "F (full model): p-value", "F statistic (proj model)",
                     "F (proj model): p-value")
    } else {
      gof.names <- c(gof.names, "F statistic", "F p-value")
    }
  }
  if (include.groups == TRUE && length(s$fe) > 0) {
    grp <- sapply(s$fe, function(x) length(levels(x)))
    grp.names <- paste0("Num. groups: ", names(grp))
    gof <- c(gof, grp)
    gof.names <- c(gof.names, grp.names)
    gof.decimal <- c(gof.decimal, rep(FALSE, length(grp)))
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

#' \code{\link{extract}} method for \code{felm} objects
#'
#' \code{\link{extract}} method for \code{felm} objects created by the
#' \code{\link[lfe]{felm}} function in the \pkg{lfe} package.
#'
#' @param model A statistical model object.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.rsquared Report R^2 in the GOF block?
#' @param include.adjrs Report adjusted R^2 in the GOF block?
#' @param include.fstatistic Report the F-statistic in the GOF block?
#' @param include.proj.stats Include statistics for projected model in the GOF block?
#' @param include.groups Report the number of groups?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract felm
#' @aliases extract.felm
#' @author Christoph Riedl, Claudia Zucca, Oliver Reiter, Philip Leifeld
#' @export
setMethod("extract", signature = className("felm", "lfe"),
          definition = extract.felm)


# -- extract.forecast (forecast) -----------------------------------------------

#' @noRd
extract.forecast <- function (model, ...) {
  model <- model$model
  return(extract(model))
}

#' \code{\link{extract}} method for \code{forecast} objects
#'
#' \code{\link{extract}} method for \code{forecast} objects created by the
#' \code{\link[forecast]{forecast}} and \code{\link[forecast:ses]{holt}} functions
#' in the \pkg{forecast} package.
#'
#' @param model A statistical model object.
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract forecast
#' @aliases extract.forecast
#' @export
setMethod("extract", signature = className("forecast", "forecast"),
          definition = extract.forecast)


# -- extract.feglm (alpaca) ----------------------------------------------------

#' @noRd
extract.feglm <- function(model, include.deviance = TRUE, include.nobs = TRUE,
                          include.groups = TRUE, ...) {
  s <- summary(model, ...)
  coefficient.names <- rownames(s$cm)
  co <- s$cm[, 1]
  se <- s$cm[, 2]
  pval <- s$cm[, 4]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()

  if (include.deviance == TRUE) {
    dev <- s$deviance
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    n <- s$nobs["nobs"]
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.groups == TRUE) {
    grp <- s$lvls.k
    grp.names <- paste0("Num. groups: ", names(grp))
    gof <- c(gof, grp)
    gof.names <- c(gof.names, grp.names)
    gof.decimal <- c(gof.decimal, rep(FALSE, length(grp)))
  }

  tr <- createTexreg(
    coef.names = coefficient.names,
    coef = co,
    se = se,
    pvalues = pval,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
  return(tr)
}

#' \code{\link{extract}} method for \code{feglm} objects
#'
#' \code{\link{extract}} method for \code{feglm} objects created by the
#' \code{\link[alpaca]{feglm}} function in the \pkg{alpaca} package.
#'
#' @param model A statistical model object.
#' @param include.deviance Report the deviance?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.groups Report the number of groups?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract feglm
#' @aliases extract.feglm
#' @author Christoph Riedl, Oliver Reiter, Philip Leifeld
#' @export
setMethod("extract", signature = className("feglm", "alpaca"),
          definition = extract.feglm)


# -- extract.fixest (fixest) ---------------------------------------------------

#' @noRd
extract.fixest <- function(model,
                           include.nobs = TRUE,       # both
                           include.groups = TRUE,     # both
                           include.rsquared = TRUE,   # OLS only
                           include.adjrs = TRUE,      # OLS only
                           include.proj.stats = TRUE, # OLS only
                           include.deviance = TRUE,   # GLM/MLE only
                           include.loglik = TRUE,     # GLM/MLE only
                           include.pseudors = TRUE,   # GLM/MLE only
                           ...) {
  # coefficient block
  s <- fixest::coeftable(model, ...)
  nam <- rownames(s)
  co <- s[, 1]
  se <- s[, 2]
  pval <- s[, 4]

  # GOF block: shared OLS and GLM/MLE statistics
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (isTRUE(include.nobs)) {
    gof <- c(gof, model$nobs)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (isTRUE(include.groups) && any(model$fixef_sizes > 0)) {
    grp <- model$fixef_sizes
    grp.names <- paste0("Num. groups: ", names(grp))
    gof <- c(gof, grp)
    gof.names <- c(gof.names, grp.names)
    gof.decimal <- c(gof.decimal, rep(FALSE, length(grp)))
  }

  # GOF block: OLS-specific statistics
  if (model$method == "feols" && isTRUE(include.rsquared)) {
    gof <- c(gof, suppressWarnings(fixest::r2(model, "r2")))
    gof.decimal <- c(gof.decimal, TRUE)
    if (isTRUE(include.proj.stats)) {
      tryCatch({ # wrap in tryCatch because it does not work in some versions
        gof <- c(gof, suppressWarnings(fixest::r2(model, "wr2")))
        gof.decimal <- c(gof.decimal, TRUE)
        gof.names <- c(gof.names, "R$^2$ (full model)", "R$^2$ (proj model)")
      }, error = function(cond) {
        gof.names <- c(gof.names, "R$^2$")
      }, warning = function(cond) {
        # do nothing in case of a warning
      })
    } else {
      gof.names <- c(gof.names, "R$^2$")
    }
  }
  if (model$method == "feols" && isTRUE(include.adjrs)) {
    gof <- c(gof, suppressWarnings(fixest::r2(model, "ar2")))
    gof.decimal <- c(gof.decimal, TRUE)
    if (isTRUE(include.proj.stats)) {
      tryCatch({ # wrap in tryCatch because it does not work in some versions
        gof <- c(gof, suppressWarnings(fixest::r2(model, "war2")))
        gof.decimal <- c(gof.decimal, TRUE)
        gof.names <- c(gof.names,
                       "Adj. R$^2$ (full model)",
                       "Adj. R$^2$ (proj model)")
      }, error = function(cond) {
        gof.names <- c(gof.names, "Adj. R$^2$")
      }, warning = function(cond) {
        # do nothing in case of a warning
      })
    } else {
      gof.names <- c(gof.names, "Adj. R$^2$")
    }
  }

  # GOF block: GLM/MLE-specific statistics
  if (model$method != "feols" && isTRUE(include.deviance)) {
    gof <- c(gof, model$deviance)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (model$method != "feols" && isTRUE(include.loglik)) {
    gof <- c(gof, model$loglik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (model$method != "feols" && isTRUE(include.pseudors)) {
    gof <- c(gof, model$pseudo_r2)
    gof.names <- c(gof.names, "Pseudo R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
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

#' \code{\link{extract}} method for \code{fixest} objects
#'
#' \code{\link{extract}} method for \code{fixest} objects created by the
#' model fitting functions in the \pkg{fixest} package. The method can deal with
#' OLS (fitted by \code{\link[fixest]{feols}}) and GLM/MLE models (fitted by
#' \code{\link[fixest]{feglm}} and other functions).
#'
#' @param model A statistical model object.
#' @param include.nobs Report the number of observations?
#' @param include.groups Report the number of groups?
#' @param include.rsquared Report R^2? (OLS only)
#' @param include.adjrs Report adjusted R^2? (OLS only)
#' @param include.proj.stats Include statistics for projected model? (OLS only)
#' @param include.deviance Report the deviance? (GLM/MLE only)
#' @param include.loglik Report the log likelihood? (GLM/MLE only)
#' @param include.pseudors Report Pseudo-R^2? (GLM/MLE only)
#' @param ... Custom parameters, which are handed over to the
#'  \code{\link[fixest]{coeftable}} method for the \code{fixest} object.
#'
#' @method extract fixest
#' @aliases extract.fixest
#' @author Christopher Poliquin, Philip Leifeld
#' @export
setMethod("extract", signature = className("fixest", "fixest"),
          definition = extract.fixest)


# -- extract.gam (mgcv) --------------------------------------------------------

#' @noRd
extract.gam <- function(model,
                        include.smooth = TRUE,
                        include.aic = TRUE,
                        include.bic = TRUE,
                        include.loglik = TRUE,
                        include.deviance = TRUE,
                        include.dev.expl = TRUE,
                        include.dispersion = TRUE,
                        include.rsquared = TRUE,
                        include.gcv = TRUE,
                        include.nobs = TRUE,
                        include.nsmooth = TRUE,
                        ...) {

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
    gof.names <- c(gof.names, "Log Likelihood")
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
    gof.names <- c(gof.names, "Deviance explained")
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
    gof.names <- c(gof.names, "GCV score")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    n <- nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.nsmooth == TRUE) {
    m <- s$m
    gof <- c(gof, m)
    gof.names <- c(gof.names, "Num. smooth terms")
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

#' \code{\link{extract}} method for \code{gam} objects
#'
#' \code{\link{extract}} method for \code{gam} objects created by the
#' \code{\link[mgcv]{gam}} function in the \pkg{mgcv} package.
#'
#' @param model A statistical model object.
#' @param include.smooth Report the smooth terms of a GAM? If they are
#'   reported, the EDF value is reported as the coefficient, and DF is included
#'   in parentheses (not standard errors because a chi-square test is used for
#'   the smooth terms).
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.deviance Report the deviance?
#' @param include.dev.expl Report the deviance explained?
#' @param include.dispersion Report the dispersion parameter?
#' @param include.rsquared Report R^2 in the GOF block?
#' @param include.gcv Report the GCV score?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.nsmooth Report the number of smooth terms?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract gam
#' @aliases extract.gam
#' @export
setMethod("extract", signature = className("gam", "mgcv"),
          definition = extract.gam)


# -- extract.bam (mgcv) --------------------------------------------------------

#' @noRd
extract.bam <- extract.gam

#' \code{\link{extract}} method for \code{bam} objects
#'
#' \code{\link{extract}} method for \code{bam} objects created by the
#' \code{\link[mgcv]{bam}} function in the \pkg{mgcv} package.
#'
#' @inheritParams extract,gam-method
#'
#' @method extract bam
#' @aliases extract.bam
#' @export
setMethod("extract", signature = className("bam", "mgcv"),
          definition = extract.bam)

# -- extract.gamlss (gamlss) ---------------------------------------------------

#' @noRd
extract.gamlss <- function(model,
                           robust = FALSE,
                           include.nobs = TRUE,
                           include.nagelkerke = TRUE,
                           include.gaic = TRUE,
                           ...) {

  # VCOV extraction; create coefficient block
  covmat <- suppressWarnings(stats::vcov(model, type = "all", robust = robust,
                                         ...))
  cf <- covmat$coef  # coefficients
  namesOfPars <- names(cf)  # names of coefficients
  se <- covmat$se  # standard errors
  tvalue <- cf / se
  pvalue <-  2 * stats::pt(-abs(tvalue), model$df.res)  # p-values

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
    n <- stats::nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{gamlss} objects
#'
#' \code{\link{extract}} method for \code{gamlss} objects created by the
#' \code{\link[gamlss]{gamlss}} function in the \pkg{gamlss} package.
#'
#' @param model A statistical model object.
#' @param robust If TRUE computes robust standard errors in the
#'   variance-covariance matrix.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.nagelkerke Report Nagelkerke R^2 in the GOF block?
#' @param include.gaic Report Generalized Akaike's Information Criterion (AIC)
#'   in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{vcov} method for the object.
#'
#' @method extract gamlss
#' @aliases extract.gamlss
#' @importFrom stats nobs pt vcov
#' @export
setMethod("extract", signature = className("gamlss", "gamlss"),
          definition = extract.gamlss)


# -- extract.gamlssZadj (gamlss.inf) -------------------------------------------

#' @noRd
extract.gamlssZadj <- function(model,
                               type = c("qr", "vcov"),
                               include.nobs = TRUE,
                               include.gaic = TRUE,
                               ...) {

  type <- match.arg(type)
  # VCOV extraction; create coefficient block
  if (type == "vcov") {
    covmat <- suppressWarnings(stats::vcov(model, type = "all", ...))
    cf <- covmat$coef  # coefficients
    namesOfPars <- names(cf)  # names of coefficients
    se <- covmat$se  # standard errors
  }
  if (type == "qr") {
    invisible(utils::capture.output(covmat <- summary(model, type = "qr")))
    cf <- covmat[, 1]
    namesOfPars <- row.names(covmat)
    se <- covmat[, 2]
  }
  tvalue <- cf / se
  pvalue <-  2 * pt(-abs(tvalue), model$df.res)  # p-values
  # add the parameter names to coefficients
  possiblePars <- c("$\\mu$", "$\\sigma$", "$\\nu$", "$\\tau$",
                    "$\\mu$ (Zero model)")
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
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
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

#' \code{\link{extract}} method for \code{gamlssZadj} objects
#'
#' \code{\link{extract}} method for \code{gamlssZadj} objects created by the
#' \code{\link[gamlss.inf]{gamlssZadj}} function in the \pkg{gamlss.inf}
#' package.
#'
#' @param model A statistical model object.
#' @param type The type.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.gaic Report Generalized Akaike's Information Criterion (AIC)
#'   in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{vcov} method for the object.
#'
#' @method extract gamlssZadj
#' @aliases extract.gamlssZadj
#' @author Ricardo Graiff Garcia, Philip Leifeld
#' @importFrom stats nobs pt vcov
#' @export
setMethod("extract", signature = className("gamlssZadj", "gamlss.inf"),
          definition = extract.gamlssZadj)


# -- extract.gee (gee) ---------------------------------------------------------

#' @noRd
extract.gee <- function(model,
                        robust = TRUE,
                        include.scale = TRUE,
                        include.nobs = TRUE,
                        ...) {
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
  pval <- 2 * stats::pnorm(abs(zval), lower.tail = FALSE)

  n <- stats::nobs(model)
  disp <- s$scale

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.scale == TRUE) {
    gof <- c(gof, disp)
    gof.names <- c(gof.names, "Scale")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{gee} objects
#'
#' \code{\link{extract}} method for \code{gee} objects created by the
#' \code{\link[gee]{gee}} function in the \pkg{gee} package.
#'
#' @param model A statistical model object.
#' @param robust If TRUE computes robust standard errors in the
#'   variance-covariance matrix.
#' @param include.scale Report the dispersion or scale parameter?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract gee
#' @aliases extract.gee
#' @importFrom stats pnorm
#' @export
setMethod("extract", signature = className("gee", "gee"),
          definition = extract.gee)


# -- extract.gel (gmm) ---------------------------------------------------------

#' @noRd
extract.gel <- function (model,
                         include.obj.fcn = TRUE,
                         include.overidentification = FALSE,
                         include.nobs = TRUE,
                         overIdentTest = c("LR", "LM", "J "),
                         ...) {

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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{gel} objects
#'
#' \code{\link{extract}} method for \code{gel} objects created by the
#' \code{\link[gmm]{gel}} function in the \pkg{gmm} package.
#'
#' @param model A statistical model object.
#' @param include.obj.fcn Report the value of the objective function
#'   (= criterion function)? More precisely, this returns
#'   \code{E(g)var(g)^{-1}E(g)}.
#' @param include.overidentification Report the J-test for overidentification?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param overIdentTest Which test statistics should be included in an
#'   overidentification test?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract gel
#' @aliases extract.gel
#' @export
setMethod("extract", signature = className("gel", "gmm"),
          definition = extract.gel)


# -- extract.geeglm (geepack) --------------------------------------------------

#' @noRd
extract.geeglm <- function(model,
                           include.scale = TRUE,
                           include.correlation = TRUE,
                           include.nobs = TRUE,
                           ...) {
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
    gof.names = c(gof.names, "Num. obs.", "Num. clust.")
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

#' \code{\link{extract}} method for \code{geeglm} objects
#'
#' \code{\link{extract}} method for \code{geeglm} objects created by the
#' \code{geeglm} function in the \pkg{geepack} package.
#'
#' @param model A statistical model object.
#' @param include.scale Report the dispersion or scale parameter?
#' @param include.correlation Report the correlation parameter alpha and its
#'   standard error?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract geeglm
#' @aliases extract.geeglm
#' @export
setMethod("extract", signature = className("geeglm", "geepack"),
          definition = extract.geeglm)


# -- extract.glm (stats) -------------------------------------------------------

#' @noRd
extract.glm <- function(model, include.aic = TRUE, include.bic = TRUE,
                        include.loglik = TRUE, include.deviance = TRUE,
                        include.nobs = TRUE, ...) {
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{glm} objects
#'
#' \code{\link{extract}} method for \code{glm} objects created by the
#' \code{\link[stats]{glm}} function in the \pkg{stats} package.
#'
#' @param model A statistical model object.
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.deviance Report the deviance?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract glm
#' @aliases extract.glm
#' @importFrom stats nobs AIC BIC logLik deviance
#' @export
setMethod("extract", signature = className("glm", "stats"),
          definition = extract.glm)


# -- extract.brglm (brglm) -----------------------------------------------------

#' @noRd
extract.brglm <- extract.glm

#' \code{\link{extract}} method for \code{brglm} objects
#'
#' \code{\link{extract}} method for \code{brglm} objects created by the
#' \code{\link[brglm]{brglm}} function in the \pkg{brglm} package.
#'
#' @inheritParams extract,glm-method
#'
#' @method extract brglm
#' @aliases extract.brglm
#' @importFrom stats pnorm
#' @export
setMethod("extract", signature = className("brglm", "brglm"),
          definition = extract.glm)


# -- extract.glm.cluster (miceadds) --------------------------------------------

#' @noRd
extract.glm.cluster <- function (model,
                                 include.aic = TRUE,
                                 include.bic = TRUE,
                                 include.loglik = TRUE,
                                 include.deviance = TRUE,
                                 include.nobs = TRUE,
                                 ...) {
  glm_model <- model$glm_res
  coefficient.names <- names(glm_model$coefficients)
  coefficients <- glm_model$coefficients
  standard.errors <- sqrt(diag(model$vcov))
  significance <- 2 * pnorm(-abs(coefficients / standard.errors))

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    aic <- AIC(glm_model)
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    bic <- BIC(glm_model)
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    lik <- logLik(glm_model)[1]
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    dev <- deviance(glm_model)
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    n <- nobs(glm_model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  tr <- createTexreg(coef.names = coefficient.names,
                     coef = coefficients,
                     se = standard.errors,
                     pvalues = significance,
                     gof.names = gof.names,
                     gof = gof,
                     gof.decimal = gof.decimal)
  return(tr)
}

#' \code{\link{extract}} method for \code{glm.cluster} objects
#'
#' \code{\link{extract}} method for \code{glm.cluster} objects created by the
#' \code{\link[miceadds:lm.cluster]{glm.cluster}} function in the \pkg{miceadds} package.
#'
#' @param model A statistical model object.
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.deviance Report the deviance?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract glm.cluster
#' @aliases extract.glm.cluster
#' @author Alexander Staudt, Philip Leifeld
#' @importFrom stats nobs AIC BIC logLik deviance
#' @export
setMethod("extract", signature = className("glm.cluster", "miceadds"),
          definition = extract.glm.cluster)


# -- extract.glmmadb (glmmADB) -------------------------------------------------

#' @noRd
extract.glmmadmb <- function(model,
                             include.variance = TRUE,
                             include.dispersion = TRUE,
                             include.zero = TRUE,
                             include.aic = TRUE,
                             include.bic = TRUE,
                             include.loglik = TRUE,
                             include.nobs = TRUE,
                             include.groups = TRUE,
                             ...) {

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
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.groups == TRUE && !is.null(model$q)) {
    groups <- model$q
    for (i in 1:length(groups)) {
      gof <- c(gof, groups[i])
      gof.names <- c(gof.names, paste("Num. groups:", names(groups)[i]))
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

#' \code{\link{extract}} method for \code{glmmadmb} objects
#'
#' \code{\link{extract}} method for \code{glmmadmb} objects created by the
#' \code{glmmadmb} function in the \pkg{glmmADMB} package.
#'
#' @param model A statistical model object.
#' @param include.variance Report group variances?
#' @param include.dispersion Report the dispersion parameter?
#' @param include.zero Should the binary part of a zero-inflated regression
#'   model or hurdle model be included in the coefficients block (after the
#'   count model)?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.groups Report the number of groups?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract glmmadmb
#' @aliases extract.glmmadmb
#' @importFrom stats pnorm
#' @export
setMethod("extract", signature = className("glmmadmb", "glmmADMB"),
          definition = extract.glmmadmb)



# -- extract.glmmTMB (glmmTMB) -------------------------------------------------

#' @noRd
extract.glmmTMB <- function(model, beside = FALSE, include.count = TRUE,
                            include.zero = TRUE, include.aic = TRUE,
                            include.groups = TRUE, include.variance = TRUE,
                            include.loglik = TRUE, include.nobs = TRUE, ...) {
  s <- summary(model, ...)
  if (model$modelInfo$allForm$ziformula == '~0') {
    include.zero <- FALSE
  }
  if (length(model$modelInfo$reTrms$cond$flist) == 0) {
    include.groups <- FALSE
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
    lik <- logLik(model)[1]
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    n <- nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.groups == TRUE) {
    grps <- sapply(model$modelInfo$reTrms$cond$flist,
                   function(x) length(levels(x)))
    grp.names <- names(grps)
    grp.names <- paste("Num. groups:", grp.names)
    gof <- c(gof, grps)
    gof.names <- c(gof.names, grp.names)
    gof.decimal <- c(gof.decimal, rep(FALSE, length(grps)))
  }
  if (include.variance == TRUE && !is.na(s$ngrps[1]) && (length(s$ngrps$cond) > 0 || length(s$ngrps$zi) > 0)) {
    vc <- glmmTMB::VarCorr(model)
    vc <- as.data.frame(rapply(vc, function(x) attr(x, "stddev")))^2
    rownames(vc) <- gsub("\\.", " ", rownames(vc))
    if (include.zero == TRUE) {
      for (i in grep("cond", rownames(vc), value = TRUE)) {
        gof.names <- c(gof.names,
                       paste("Var (count model):", sub("cond ", "", i)))
      }
      for (i in grep("zi", rownames(vc), value = TRUE)) {
        gof.names <- c(gof.names, paste("Var (zero model):", sub("zi ", "", i)))
      }
    } else {
      for (i in grep("cond", rownames(vc), value = TRUE)) {
        gof.names <- c(gof.names, paste("Var:", sub("cond ", "", i)))
      }
    }
    gof <- c(gof, vc[, 1])
    gof.decimal <- c(gof.decimal, rep(TRUE, nrow(vc)))
  }
  count <- coef(s)$cond
  zero <- coef(s)$zi
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
      stop(paste("Either the 'include.count' or the 'include.zero' argument",
                 "must be TRUE."))
    }
    names <- rownames(coef.block)
    co <- coef.block[, 1]
    se <- coef.block[, 2]
    pval <- coef.block[, 4]
    tr <- createTexreg(coef.names = names,
                       coef = co, se = se,
                       pvalues = pval,
                       gof.names = gof.names,
                       gof = gof,
                       gof.decimal = gof.decimal)
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
      tr <- createTexreg(coef.names = c.names,
                         coef = c.co,
                         se = c.se,
                         pvalues = c.pval,
                         gof.names = gof.names,
                         gof = gof,
                         gof.decimal = gof.decimal,
                         model.name = "Count model")
      trList[[length(trList) + 1]] <- tr
    }
    if (include.zero == TRUE) {
      tr <- createTexreg(coef.names = z.names,
                         coef = z.co,
                         se = z.se,
                         pvalues = z.pval,
                         gof.names = gof.names,
                         gof = gof,
                         gof.decimal = gof.decimal,
                         model.name = "Zero model")
      trList[[length(trList) + 1]] <- tr
    }
    if (length(trList) == 0) {
      stop(paste("Either the 'include.count' or the 'include.zero' argument",
                 "must be TRUE."))
    }
    return(trList)
  }
}

#' \code{\link{extract}} method for \code{glmmTMB} objects
#'
#' \code{\link{extract}} method for \code{glmmTMB} objects created by the
#' \code{\link[glmmTMB]{glmmTMB}} function in the \pkg{glmmTMB} package.
#'
#' @param model A statistical model object.
#' @param beside Arrange the model terms below each other or beside each other?
#'   The binary model parameters and the count parameters can be displayed in
#'   two separate columns of the table.
#' @param include.count Report the count parameters in the coefficients block
#'   (before the binary part for the zeros)?
#' @param include.zero Should the binary part of the model be included in the
#'   coefficients block (after the count parameters)?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.groups Report the number of groups?
#' @param include.variance Report group variances?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract glmmTMB
#' @aliases extract.glmmTMB
#' @author Ricardo Graiff Garcia, Philip Leifeld
#' @importFrom stats AIC logLik nobs
#' @export
setMethod("extract", signature = className("glmmTMB", "glmmTMB"),
          definition = extract.glmmTMB)


# -- extract.gls (nlme) --------------------------------------------------------

#' @noRd
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{gls} objects
#'
#' \code{\link{extract}} method for \code{gls} objects created by the
#' \code{\link[nlme]{gls}} function in the \pkg{nlme} package.
#'
#' @param model A statistical model object.
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract gls
#' @aliases extract.gls
#' @export
setMethod("extract", signature = className("gls", "nlme"),
          definition = extract.gls)


# -- extract.gmm (gmm) -------------------------------------------------------

#' @noRd
extract.gmm <- function(model,
                        include.obj.fcn = TRUE,
                        include.overidentification = FALSE,
                        include.nobs = TRUE,
                        ...) {

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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{gmm} objects
#'
#' \code{\link{extract}} method for \code{gmm} objects created by the
#' \code{\link[gmm]{gmm}} function in the \pkg{gmm} package.
#'
#' @inheritParams extract,gel-method
#'
#' @method extract gmm
#' @aliases extract.gmm
#' @export
setMethod("extract", signature = className("gmm", "gmm"),
          definition = extract.gmm)


# -- extract.gnls (nlme) -------------------------------------------------------

#' @noRd
extract.gnls <- extract.gls

#' \code{\link{extract}} method for \code{gnls} objects
#'
#' \code{\link{extract}} method for \code{gnls} objects created by the
#' \code{\link[nlme]{gnls}} function in the \pkg{nlme} package.
#'
#' @inheritParams extract,gls-method
#'
#' @method extract gnls
#' @aliases extract.gnls
#' @export
setMethod("extract", signature = className("gnls", "nlme"),
          definition = extract.gnls)


# -- extract.gnm (gnm) ---------------------------------------------------------

#' @noRd
extract.gnm <- function(model, include.aic = TRUE, include.bic = TRUE,
                        include.loglik = TRUE, include.deviance = TRUE,
                        include.nobs = TRUE, include.df = FALSE,
                        include.chisq = FALSE, include.delta = FALSE, ...) {

  s <- summary(model)
  coefficients.names <- names(model$coefficients)
  co <- s$coef[, 1]
  se <- s$coef[, 2]
  pval <- s$coef[, 3]

  table <- as.data.frame(cbind(coefficients.names, co, se, pval))
  table <- table[!is.na(table$se), ]
  coefficients.names <- as.character(table$coefficients.names)
  co <- as.numeric(as.character(table$co))
  se <- as.numeric(as.character(table$se))
  pval <- as.numeric(as.character(table$pval))

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()

  if (include.df == TRUE) {
    df <- model$df.residual
    gof <- c(gof, df)
    gof.names <- c(gof.names, "df. residuals")
    gof.decimal <- c(gof.decimal, FALSE)
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
  if (include.chisq == TRUE) {
    chisq <- sum(na.omit(c(residuals(model, "pearson")^2)))
    gof <- c(gof, chisq)
    gof.names <- c(gof.names, "Pearson chi-squared")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.aic == TRUE) {
    aic <- AIC(model)[1]
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    bic <- BIC(model)[1]
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.delta == TRUE) {
    delta <- sum(na.omit(c(abs(residuals(model,"response"))))) /
      sum(na.omit(c(abs(fitted(model))))) / 2 * 100 # Dissimilarity index
    gof <- c(gof, delta)
    gof.names <- c(gof.names, "Dissim. Index")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    nobs <- length(model$y)
    gof <- c(gof, nobs)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }

  tr <- createTexreg(
    coef.names = coefficients.names,
    coef = co,
    se = se,
    pvalues = pval,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
  return(tr)
}

#' \code{\link{extract}} method for \code{gnm} objects
#'
#' \code{\link{extract}} method for \code{gnm} objects created by the
#' \code{\link[gnm]{gnm}} function in the \pkg{gnm} package.
#'
#' @param model A statistical model object.
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.deviance Report the deviance?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.df Report the degrees of freedom?
#' @param include.chisq Report the chi squared statistic?
#' @param include.delta Report the delta statistic?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract gnm
#' @aliases extract.gnm
#' @importFrom stats nobs AIC BIC residuals logLik na.omit fitted
#' @export
setMethod("extract", signature = className("gnm", "gnm"),
          definition = extract.gnm)


# -- extract.H2OBinomialModel (h20) --------------------------------------------

#' @noRd
extract.H2OBinomialModel <- function(model,
                                     standardized = FALSE,
                                     include.mse = TRUE,
                                     include.rsquared = TRUE,
                                     include.logloss = TRUE,
                                     include.meanerror = TRUE,
                                     include.auc = TRUE,
                                     include.gini = TRUE,
                                     include.deviance = TRUE,
                                     include.aic = TRUE,
                                     ...) {

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

#' \code{\link{extract}} method for \code{H2OBinomialModel} objects
#'
#' \code{\link{extract}} method for \code{H2OBinomialModel} objects created by
#' the \code{\link[h2o]{h2o.glm}} function in the \pkg{h2o} package.
#'
#' @param model A statistical model object.
#' @param standardized Report standardized coefficients instead of raw
#'   coefficients?
#' @param include.mse Report the mean squared error in the GOF block?
#' @param include.rsquared Report R^2 in the GOF block?
#' @param include.logloss Report the log loss?
#' @param include.meanerror Report the mean per-class error?
#' @param include.auc Report the area under the curve (AUC)?
#' @param include.gini Report the Gini coefficient?
#' @param include.deviance Report the deviance?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract H2OBinomialModel
#' @aliases extract.H2OBinomialModel
#' @export
setMethod("extract", signature = className("H2OBinomialModel", "h2o"),
          definition = extract.H2OBinomialModel)


# -- extract.lm (stats) --------------------------------------------------------

#' @noRd
extract.lm <- function(model, include.rsquared = TRUE, include.adjrs = TRUE,
                       include.nobs = TRUE, include.fstatistic = FALSE,
                       include.rmse = FALSE, ...) {
  s <- summary(model, ...)

  names <- rownames(s$coefficients)
  co <- s$coefficients[, 1]
  se <- s$coefficients[, 2]
  pval <- s$coefficients[, 4]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (isTRUE(include.rsquared)) {
    rs <- s$r.squared  # extract R-squared
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (isTRUE(include.adjrs)) {
    adj <- s$adj.r.squared  # extract adjusted R-squared
    gof <- c(gof, adj)
    gof.names <- c(gof.names, "Adj. R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (isTRUE(include.nobs)) {
    n <- nobs(model)  # extract number of observations
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (isTRUE(include.fstatistic)) {
    fstat <- s$fstatistic[[1]]
    gof <- c(gof, fstat)
    gof.names <- c(gof.names, "F statistic")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (isTRUE(include.rmse) && !is.null(s$sigma[[1]])) {
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

#' \code{\link{extract}} method for \code{lm} objects
#'
#' \code{\link{extract}} method for \code{lm} objects created by the
#' \code{\link[stats]{lm}} function in the \pkg{stats} package.
#'
#' @param model A statistical model object.
#' @param include.rsquared Report R^2 in the GOF block?
#' @param include.adjrs Report adjusted R^2 in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.fstatistic Report the F-statistic in the GOF block?
#' @param include.rmse Report the root mean square error (RMSE; = residual
#'   standard deviation) in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract lm
#' @aliases extract.lm
#' @importFrom stats nobs
#' @export
setMethod("extract", signature = className("lm", "stats"),
          definition = extract.lm)


# -- extract.dynlm (dynlm) -----------------------------------------------------

#' @noRd
extract.dynlm <- extract.lm

#' \code{\link{extract}} method for \code{dynlm} objects
#'
#' \code{\link{extract}} method for \code{dynlm} objects created by the
#' \code{\link[dynlm]{dynlm}} function in the \pkg{dynlm} package.
#'
#' @inheritParams extract,lm-method
#'
#' @method extract dynlm
#' @aliases extract.dynlm
#' @export
setMethod("extract", signature = className("dynlm", "dynlm"),
          definition = extract.dynlm)


# -- extract.ivreg (AER) -------------------------------------------------------

#' @noRd
extract.ivreg <- extract.lm

#' \code{\link{extract}} method for \code{ivreg} objects
#'
#' \code{\link{extract}} method for \code{ivreg} objects created by the
#' \code{\link[AER]{ivreg}} function in the \pkg{AER} package.
#'
#' @inheritParams extract,lm-method
#'
#' @method extract ivreg
#' @aliases extract.ivreg
#' @export
setMethod("extract", signature = className("ivreg", "AER"),
          definition = extract.ivreg)


# -- extract.lm.cluster (miceadds) ---------------------------------------------

#' @noRd
extract.lm.cluster <- function(model,
                               include.rsquared = TRUE,
                               include.adjrs = TRUE,
                               include.nobs = TRUE,
                               include.fstatistic = FALSE,
                               include.rmse = FALSE,
                               ...) {

  s <- summary(model$lm_res)
  lm_model <- model$lm_res
  nam <- names(lm_model$coefficients)
  co <- lm_model$coefficients
  se <- sqrt(diag(model$vcov))
  pval <- 2 * pnorm(-abs(co / se))


  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.rsquared == TRUE) {
    rs <- s$r.squared  # extract R-squared
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.adjrs == TRUE) {
    adj <- s$adj.r.squared  # extract adjusted R-squared
    gof <- c(gof, adj)
    gof.names <- c(gof.names, "Adj. R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    n <- nobs(lm_model)  # extract number of observations
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{lm.cluster} objects
#'
#' \code{\link{extract}} method for \code{lm.cluster} objects created by the
#' \code{\link[miceadds]{lm.cluster}} function in the \pkg{miceadds} package.
#'
#' @param model A statistical model object.
#' @param include.rsquared Report R^2 in the GOF block?
#' @param include.adjrs Report adjusted R^2 in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.fstatistic Report the F-statistic in the GOF block?
#' @param include.rmse Report the root mean square error (RMSE; = residual
#'   standard deviation) in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract lm.cluster
#' @aliases extract.lm.cluster
#' @author Alexander Staudt, Philip Leifeld
#' @export
setMethod("extract", signature = className("lm.cluster", "miceadds"),
          definition = extract.lm.cluster)


# -- extract.lme (nlme) --------------------------------------------------------

#' @noRd
extract.lme <- function(model, include.aic = TRUE, include.bic = TRUE,
                        include.loglik = TRUE, include.nobs = TRUE,
                        include.groups = TRUE, include.variance = FALSE, ...) {

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
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.groups == TRUE) {
    grp <- model$dims$ngrps[1:model$dims$Q]
    for(i in 1:length(grp)){
      gof <- c(gof, grp[i])
      gof.names <- c(gof.names, paste("Num. groups:", names(grp)[i]))
      gof.decimal <- c(gof.decimal, FALSE)
    }
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
        gof.names <- c(gof.names, "sigma. RE")
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

#' \code{\link{extract}} method for \code{lme} objects
#'
#' \code{\link{extract}} method for \code{lme} objects created by the
#' \code{\link[nlme]{lme}} function in the \pkg{nlme} package.
#'
#' @param model A statistical model object.
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.groups Report the number of groups?
#' @param include.variance Report group variances?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract lme
#' @aliases extract.lme
#' @export
setMethod("extract", signature = className("lme", "nlme"),
          definition = extract.lme)


# -- extract.glmmPQL (MASS) ----------------------------------------------------

#' @noRd
extract.glmmPQL <- extract.lme

#' \code{\link{extract}} method for \code{glmmPQL} objects
#'
#' \code{\link{extract}} method for \code{glmmPQL} objects created by the
#' \code{\link[MASS]{glmmPQL}} function in the \pkg{MASS} package.
#'
#' @inheritParams extract,lme-method
#'
#' @method extract glmmPQL
#' @aliases extract.glmmPQL
#' @export
setMethod("extract", signature = className("glmmPQL", "MASS"),
          definition = extract.glmmPQL)


# -- extract.lme4 (lme4) -------------------------------------------------------

#' @noRd
extract.lme4 <- function(model,
                         method = c("naive", "profile", "boot", "Wald"),
                         level = 0.95,
                         nsim = 1000,
                         include.aic = TRUE,
                         include.bic = TRUE,
                         include.dic = FALSE,
                         include.deviance = FALSE,
                         include.loglik = TRUE,
                         include.nobs = TRUE,
                         include.groups = TRUE,
                         include.variance = TRUE,
                         ...) {

  if (utils::packageVersion("lme4") < "1.0") {
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
    n <- dim(stats::model.frame(model))[1]
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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
    # vc <- lme4::VarCorr(model)
    # varcomps <- c(unlist(lapply(vc, diag)),   # random intercept variances
    #               attr(vc, "sc")^2)                     # residual variance
    # varnames <- names(varcomps)
    # varnames[length(varnames)] <- "Residual"
    # varnames <- paste("Variance:", varnames)
    # if (is.na(attr(vc, "sc"))) {
    #   varnames <- varnames[-length(varnames)]
    #   varcomps <- varcomps[-length(varcomps)]
    # }
    # gof <- c(gof, varcomps)
    # gof.names <- c(gof.names, varnames)
    # gof.decimal <- c(gof.decimal, rep(TRUE, length(varcomps)))
  }

  betas <- lme4::fixef(model, ...)
  if ("confint.merMod" %in% utils::methods("confint") && method[1] != "naive") {
    ci <- tryCatch({
      ci <- stats::confint(model, method = method[1], level = level, nsim = nsim, ...)
    },
    error = function(err) {
      method <- "naive"
      message("Confidence intervals not available for this model. Using naive p-values instead.")
    }
    )
    if (is.null(ci)) {
      method <- "naive"
    } else {
      last <- nrow(ci)
      number <- length(betas)
      first <- last - number + 1
      ci <- ci[first:last, ]
      if (is.matrix(ci)) {
        ci.l <- ci[, 1]
        ci.u <- ci[, 2]
      } else {
        ci.l <- ci[1]
        ci.u <- ci[2]
      }
    }
  } else if (method[1] != "naive") {
    method[1] <- "naive"
    message("confint.merMod method not found. Using naive p-values instead.")
  }

  if (method[1] == "naive") {
    Vcov <- tryCatch({
      Vcov <- stats::vcov(model, useScale = FALSE, ...)
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

#' \code{\link{extract}} method for \code{lme4} objects
#'
#' \code{\link{extract}} method for \code{lme4} objects created by the
#' \pkg{lme4} package.
#'
#' @param model A statistical model object.
#' @param method The method used to compute confidence intervals or p-values.
#'   The default value \code{"naive"} computes naive p-values while the other
#'   methods compute confidence intervals using the \code{confint} function. See
#'   \code{\link[lme4]{confint.merMod}}.
#' @param level Significance or confidence level (\code{1 - alpha}) for
#'   computing confidence intervals.
#' @param nsim The MCMC sample size or number of bootstrapping replications on
#'   the basis of which confidence intervals are computed (only if the
#'   \code{method} argument does not specify \code{"naive"}, which is the
#'   default behavior). Note: large values may take considerable computing time.
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.dic Report the deviance information criterion (DIC)?
#' @param include.deviance Report the deviance?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.groups Report the number of groups?
#' @param include.variance Report group variances?
#' @param ... Arguments to be passed to the \code{\link[lme4]{fixef}} function
#'   in the \pkg{lme4} package.
#'
#' @method extract lme4
#' @aliases extract.lme4
#' @importFrom stats confint model.frame vcov
#' @importFrom utils methods packageVersion
#' @export
setMethod("extract", signature = className("lme4", "lme4"),
          definition = extract.lme4)


# -- extract.merMod (lme4) -----------------------------------------------------

#' @noRd
extract.merMod <- extract.lme4

#' \code{\link{extract}} method for \code{merMod} objects
#'
#' \code{\link{extract}} method for \code{merMod} objects created by the
#' \pkg{lme4} package.
#'
#' @inheritParams extract,lme4-method
#'
#' @method extract merMod
#' @aliases extract.merMod
#' @export
setMethod("extract", signature = className("merMod", "lme4"),
          definition = extract.merMod)


# -- extract.lmerMod (lme4) ----------------------------------------------------

#' @noRd
extract.lmerMod <- extract.lme4

#' \code{\link{extract}} method for \code{lmerMod} objects
#'
#' \code{\link{extract}} method for \code{lmerMod} objects created by the
#' \code{\link[lme4]{lmer}} function in the \pkg{lme4} package.
#'
#' @inheritParams extract,lme4-method
#'
#' @method extract lmerMod
#' @aliases extract.lmerMod
#' @export
setMethod("extract", signature = className("lmerMod", "lme4"),
          definition = extract.lmerMod)


# -- extract.glmerMod (lme4) ---------------------------------------------------

#' @noRd
extract.glmerMod <- extract.lme4

#' \code{\link{extract}} method for \code{glmerMod} objects
#'
#' \code{\link{extract}} method for \code{glmerMod} objects created by the
#' \code{\link[lme4]{glmer}} function in the \pkg{lme4} package.
#'
#' @inheritParams extract,lme4-method
#'
#' @method extract glmerMod
#' @aliases extract.glmerMod
#' @export
setMethod("extract", signature = className("glmerMod", "lme4"),
          definition = extract.glmerMod)


# -- extract.nlmerMod (lme4) ---------------------------------------------------

#' @noRd
extract.nlmerMod <- extract.lme4

#' \code{\link{extract}} method for \code{nlmerMod} objects
#'
#' \code{\link{extract}} method for \code{nlmerMod} objects created by the
#' \code{\link[lme4]{nlmer}} function in the \pkg{lme4} package.
#'
#' @inheritParams extract,lme4-method
#'
#' @method extract nlmerMod
#' @aliases extract.nlmerMod
#' @export
setMethod("extract", signature = className("nlmerMod", "lme4"),
          definition = extract.nlmerMod)


# -- extract.lmrob (robustbase) ------------------------------------------------

#' @noRd
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{lmrob} objects
#'
#' \code{\link{extract}} method for \code{lmrob} objects created by the
#' \code{lmrob} function in the \pkg{robustbase} package.
#'
#' @param model A statistical model object.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract lmrob
#' @aliases extract.lmrob
#' @export
setMethod("extract", signature = className("lmrob", "robustbase"),
          definition = extract.lmrob)


# -- extract.glmrob (robustbase) -----------------------------------------------

#' @noRd
extract.glmrob <- extract.lmrob

#' \code{\link{extract}} method for \code{glmrob} objects
#'
#' \code{\link{extract}} method for \code{glmrob} objects created by the
#' \code{glmrob} function in the \pkg{robustbase} package.
#'
#' @inheritParams extract,lmrob-method
#'
#' @method extract glmrob
#' @aliases extract.glmrob
#' @export
setMethod("extract", signature = className("glmrob", "robustbase"),
          definition = extract.glmrob)


# -- extract.lmRob (robust) ----------------------------------------------------

#' @noRd
extract.lmRob <- function(model, include.rsquared = TRUE,
                          include.nobs = TRUE, include.rmse = TRUE, ...) {
  s <- summary(model, ...)

  names <- rownames(s$coefficients)
  co <- s$coefficients[, 1]
  se <- s$coefficients[, 2]
  pval <- s$coefficients[, 4]

  rs <- s$r.squared  #extract R-squared
  n <- length(model$residuals)  #extract number of observations

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.rsquared == TRUE) {
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.rmse == TRUE && !is.null(s$sigma)) {
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

#' \code{\link{extract}} method for \code{lmRob} objects
#'
#' \code{\link{extract}} method for \code{lmRob} objects created by the
#' \code{\link[robust]{lmRob}} function in the \pkg{robust} package.
#'
#' @param model A statistical model object.
#' @param include.rsquared Report R^2 in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.rmse Report the root mean square error (RMSE; = residual
#'   standard deviation) in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract lmRob
#' @aliases extract.lmRob
#' @rdname extract-lmrob-method
#' @export
setMethod("extract", signature = className("lmRob", "robust"),
          definition = extract.lmRob)


# -- extract.lnam (sna) --------------------------------------------------------

#' @noRd
extract.lnam <- function(model,
                         include.rsquared = TRUE,
                         include.adjrs = TRUE,
                         include.aic = TRUE,
                         include.bic = TRUE,
                         include.loglik = TRUE,
                         ...) {
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
    gof.names <- c(gof.names, "Adj. R$^2$")
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

#' \code{\link{extract}} method for \code{lnam} objects
#'
#' \code{\link{extract}} method for \code{lnam} objects created by the
#' \code{lnam} function in the \pkg{sna} package.
#'
#' @param model A statistical model object.
#' @param include.rsquared Report R^2 in the GOF block?
#' @param include.adjrs Report adjusted R^2 in the GOF block?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{coef} method for the object.
#'
#' @method extract lnam
#' @aliases extract.lnam
#' @importFrom stats pnorm
#' @export
setMethod("extract", signature = className("lnam", "sna"),
          definition = extract.lnam)


# -- extract.logitmfx (mfx) ----------------------------------------------------

#' @noRd
extract.logitmfx <- function(model,
                             include.nobs = TRUE,
                             include.loglik = TRUE,
                             include.deviance = TRUE,
                             include.aic = TRUE,
                             include.bic = TRUE,
                             ...) {
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{logitmfx} objects
#'
#' \code{\link{extract}} method for \code{logitmfx} objects created by the
#' \code{\link[mfx]{logitmfx}} function in the \pkg{mfx} package.
#'
#' @param model A statistical model object.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.deviance Report the deviance?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract logitmfx
#' @aliases extract.logitmfx
#' @export
setMethod("extract", signature = className("logitmfx", "mfx"),
          definition = extract.logitmfx)


# -- extract.probitmfx (mfx) ---------------------------------------------------

#' @noRd
extract.probitmfx <- extract.logitmfx

#' \code{\link{extract}} method for \code{probitmfx} objects
#'
#' \code{\link{extract}} method for \code{probitmfx} objects created by the
#' \code{\link[mfx]{probitmfx}} function in the \pkg{mfx} package.
#'
#' @inheritParams extract,logitmfx-method
#'
#' @method extract probitmfx
#' @aliases extract.probitmfx
#' @export
setMethod("extract", signature = className("probitmfx", "mfx"),
          definition = extract.probitmfx)


# -- extract.logitor (mfx) -----------------------------------------------------

#' @noRd
extract.logitor <- function(model,
                            include.nobs = TRUE,
                            include.loglik = TRUE,
                            include.deviance = TRUE,
                            include.aic = TRUE,
                            include.bic = TRUE,
                            ...) {
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{logitor} objects
#'
#' \code{\link{extract}} method for \code{logitor} objects created by the
#' \code{\link[mfx]{logitor}} function in the \pkg{mfx} package.
#'
#' @inheritParams extract,logitmfx-method
#'
#' @method extract logitor
#' @aliases extract.logitor
#' @export
setMethod("extract", signature = className("logitor", "mfx"),
          definition = extract.logitor)


# -- extract.lqmm (lqmm) -------------------------------------------------------

#' @noRd
extract.lqmm <- function(model,
                         include.aic = TRUE,
                         include.bic = TRUE,
                         include.loglik = TRUE,
                         include.nobs = TRUE,
                         include.groups = TRUE,
                         include.tau = FALSE,
                         use.ci = FALSE,
                         beside = TRUE,
                         ...) {

  s <- summary(model, ...)

  tau <- model$tau
  if (length(tau) == 1 && !"list" %in% class(s$tTable)) {
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
        gof.names <- c(gof.names, "Num. obs.")
        gof.decimal <- c(gof.decimal, FALSE)
      }
      if (include.groups == TRUE) {
        gof <- c(gof, model$ngroups)
        gof.names <- c(gof.names, "Num. groups")
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
      gof.names <- c(gof.names, "Num. obs.")
      gof.decimal <- c(gof.decimal, FALSE)
    }
    if (include.groups == TRUE) {
      gof <- c(gof, model$ngroups)
      gof.names <- c(gof.names, "Num. groups")
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

#' \code{\link{extract}} method for \code{lqmm} objects
#'
#' \code{\link{extract}} method for \code{lqmm} objects created by the
#' \code{\link[lqmm]{lqmm}} function in the \pkg{lqmm} package.
#'
#' @param model A statistical model object.
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.groups Report the number of groups?
#' @param include.tau Report tau?
#' @param use.ci Report confidence intervals in the GOF block?
#' @param beside Arrange the model terms below each other or beside each other?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract lqmm
#' @aliases extract.lqmm
#' @importFrom stats AIC BIC
#' @export
setMethod("extract", signature = className("lqmm", "lqmm"),
          definition = extract.lqmm)


# -- extract.lrm (rms) ---------------------------------------------------------

#' @noRd
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{lrm} objects
#'
#' \code{\link{extract}} method for \code{lrm} objects created by the
#' \code{\link[rms]{lrm}} function in the \pkg{rms} package.
#'
#' @param model A statistical model object.
#' @param include.pseudors Report pseudo R^2 in the GOF block?
#' @param include.lr Report likelihood ratio test?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract lrm
#' @aliases extract.lrm
#' @author Fabrice Le Lec
#' @importFrom stats pnorm
#' @export
setMethod("extract", signature = className("lrm", "rms"),
          definition = extract.lrm)


# -- extract.maxLik (maxLik) ---------------------------------------------------

#' @noRd
extract.maxLik <- function(model,
                           include.loglik = TRUE,
                           include.aic = TRUE,
                           ...) {

  s <- summary(model, ...)
  coefs <- s$estimate[, 1]
  se <- s$estimate[, 2]
  pval <- s$estimate[, 4]
  coefnames <- rownames(s$estimate)
  if (is.null(coefnames)) {
    coefnames <- character(length(coefs))
  }

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (isTRUE(include.loglik)) {
    ll <- stats::logLik(model)[1]
    if (!is.null(ll) && !is.na(ll)) {
      gof <- c(gof, ll)
      gof.names <- c(gof.names, "Log Likelihood")
      gof.decimal <- c(gof.decimal, TRUE)
    }
  }
  if (isTRUE(include.aic)) {
    aic <- stats::AIC(model)[1]
    if (!is.null(aic) && !is.na(aic)) {
      gof <- c(gof, aic)
      gof.names <- c(gof.names, "AIC")
      gof.decimal <- c(gof.decimal, TRUE)
    }
  }

  createTexreg(coef.names = coefnames,
               coef = coefs,
               se = se,
               pvalues = pval,
               gof.names = gof.names,
               gof = gof,
               gof.decimal = gof.decimal)
}

#' \code{\link{extract}} method for \code{maxLik} objects
#'
#' \code{\link{extract}} method for \code{maxLik} objects created by the
#' \code{\link[maxLik]{maxLik}} function in the \pkg{maxLik} package.
#'
#' @param model A statistical model object.
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.aic Report the AIC in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract maxLik
#' @aliases extract.maxLik
#' @importFrom stats AIC logLik
#' @export
setMethod("extract", signature = className("maxLik", "maxLik"),
          definition = extract.maxLik)


# -- extract.mhurdle (mhurdle) -------------------------------------------------

#' @noRd
extract.mhurdle <- function (model,
                             include.nobs = TRUE,
                             include.loglik = TRUE,
                             ...) {

  s <- summary(model, ...)
  names <- rownames(s$coefficients)
  class(names) <- "character"
  co <- s$coefficients[, 1]
  se <- s$coefficients[, 2]
  pval <- s$coefficients[, 4]
  class(co) <- class(se) <- class(pval) <- "numeric"

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.loglik == TRUE) {
    gof <- c(gof, as.numeric(s$naive$logLik))
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, length(s$model[, 1]))
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }

  tr <- createTexreg(coef.names = names,
                     coef = co,
                     se = se,
                     pvalues = pval,
                     gof.names = gof.names,
                     gof = gof,
                     gof.decimal = gof.decimal)
  return(tr)
}

#' \code{\link{extract}} method for \code{mhurdle} objects
#'
#' \code{\link{extract}} method for \code{mhurdle} objects created by the
#' \code{\link[mhurdle]{mhurdle}} function in the \pkg{mhurdle} package.
#'
#' @param model A statistical model object.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract mhurdle
#' @aliases extract.mhurdle
#' @export
setMethod("extract", signature = className("mhurdle", "mhurdle"),
          definition = extract.mhurdle)


# -- extract.mlogit (mlogit) ---------------------------------------------------

#' @noRd
extract.mlogit <- function(model,
                           include.aic = TRUE,
                           include.loglik = TRUE,
                           include.nobs = TRUE,
                           include.groups = TRUE,
                           include.order = FALSE,
                           include.iterations = FALSE,
                           beside = FALSE,
                           ...) {
  s <- summary(model, ...)

  if (include.order == TRUE) {
    s$CoefTable <- s$CoefTable[order(rownames(s$CoefTable)),]
  }

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
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, nrow(s$residuals))
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.groups == TRUE) {
    K <- ncol(s$residuals)
    gof <- c(gof, K)
    gof.names <- c(gof.names, "K")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.iterations == TRUE) {
    iter <- s$est.stat$nb.iter
    gradNorm <- s$est.stat$eps[1, 1]
    gof <- c(gof, iter, gradNorm)
    gof.names <- c(gof.names, "Iterations", "Gradient 2-norm")
    gof.decimal <- c(gof.decimal, c(FALSE, TRUE))
  }

  # check for choice-specific covariates and fix beside argument if necessary
  models <- attributes(model$freq)$dimnames[[1]][-1]
  rows <- which(grepl(paste0(":", models[1], "$"), rn))
  if (isTRUE(beside) &&
      (ncol(s$residuals) - 1) * length(rows) < nrow(s$CoefTable)) {
    beside <- FALSE
    warning(paste0("The mlogit model has choice-specific covariates; 'beside'",
                   " argument will be set to FALSE."))
  }

  if (isTRUE(beside)) {
    models <- attributes(model$freq)$dimnames[[1]][-1]
    trlist <- list()
    for (i in 1:length(models)) {
      rows <- which(grepl(paste0(":", models[i], "$"), rn))
      coeftable <- s$CoefTable[rows, ]
      cn <- rn[rows]
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
  } else {
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
}

#' \code{\link{extract}} method for \code{mlogit} objects
#'
#' \code{\link{extract}} method for \code{mlogit} objects created by the
#' \code{\link[mlogit]{mlogit}} function in the \pkg{mlogit} package.
#'
#' @param model A statistical model object.
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.groups Report the number of groups?
#' @param include.order Report coefficient names in alphabetical order?
#' @param include.iterations Report the number of iterations?
#' @param beside Arrange the model terms below each other or beside each other?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract mlogit
#' @aliases extract.mlogit
#' @importFrom stats AIC logLik
#' @export
setMethod("extract", signature = className("mlogit", "mlogit"),
          definition = extract.mlogit)


# -- extract.model.selection (MuMIn) -------------------------------------------

#' @noRd
extract.model.selection <- function(model,
                                    include.loglik = TRUE,
                                    include.aicc = TRUE,
                                    include.delta = TRUE,
                                    include.weight = TRUE,
                                    include.nobs = TRUE,
                                    ...) {

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
                             nobs = "Num. obs.")[include])

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

#' \code{\link{extract}} method for \code{model.selection} objects
#'
#' \code{\link{extract}} method for \code{model.selection} objects created by
#' the \code{\link[MuMIn]{model.sel}} and \code{\link[MuMIn]{dredge}} functions
#' in the \pkg{MuMIn} package.
#'
#' @param model A statistical model object.
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.aicc Report AICC in the GOF block?
#' @param include.delta Report the delta statistic?
#' @param include.weight Report Akaike weights?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract model.selection
#' @aliases extract.model.selection
#' @export
setMethod("extract", signature = className("model.selection", "MuMIn"),
          definition = extract.model.selection)


# -- extract.multinom (nnet) ---------------------------------------------------

# extension for multinom objects (nnet package)
extract.multinom <- function(model,
                             include.pvalues = TRUE,
                             include.aic = TRUE,
                             include.bic = TRUE,
                             include.loglik = TRUE,
                             include.deviance = TRUE,
                             include.nobs = TRUE,
                             include.groups = TRUE,
                             levels = model$lev,
                             beside = FALSE,
                             ...) {

  s <- summary(model, ...)

  coefnames <- model$coefnames
  co <- s$coefficients
  se <- s$standard.errors

  if (!"matrix" %in% class(co)) {
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
    n <- nrow(s$fitted.values)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.groups == TRUE) {
    K <- ncol(model$residuals)
    if (K == 1) {
      K <- 2
    }
    gof <- c(gof, K)
    gof.names <- c(gof.names, "K")
    gof.decimal <- c(gof.decimal, FALSE)
  }

  if (is.null(rownames(co))) {
    if (include.pvalues == TRUE) {
      zval <- co[1, ] / se[1, ]
      pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    } else {
      pval <- numeric(0)
    }
    tr <- createTexreg(
      coef.names = coefnames,
      coef = co[1, ],
      se = se[1, ],
      pvalues = pval,
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal
    )
    return(tr)
  } else if (beside == TRUE) {
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
      } else {
        pval <- numeric(0)
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

#' \code{\link{extract}} method for \code{multinom} objects
#'
#' \code{\link{extract}} method for \code{multinom} objects created by the
#' \code{\link[nnet]{multinom}} function in the \pkg{nnet} package.
#'
#' @param model A statistical model object.
#' @param include.pvalues Report p-values?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.deviance Report the deviance?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.groups Report the number of groups?
#' @param levels The names of the levels of a multinomial model that should be
#'   included in the table. Should be provided as a vector of character strings.
#' @param beside Arrange the model terms below each other or beside each other?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract multinom
#' @aliases extract.multinom
#' @importFrom stats AIC BIC logLik deviance pnorm
#' @export
setMethod("extract", signature = className("multinom", "nnet"),
          definition = extract.multinom)


# -- extract.negbin (MASS) -----------------------------------------------------

#' @noRd
extract.negbin <- extract.glm

#' \code{\link{extract}} method for \code{negbin} objects
#'
#' \code{\link{extract}} method for \code{negbin} objects created by the
#' \code{\link[MASS]{glm.nb}} function in the \pkg{MASS} package.
#'
#' @inheritParams extract,glm-method
#'
#' @method extract negbin
#' @aliases extract.negbin
#' @export
setMethod("extract", signature = className("negbin", "MASS"),
          definition = extract.negbin)


# -- extract.negbinirr (mfx) ---------------------------------------------------

#' @noRd
extract.negbinirr <- function(model,
                              include.nobs = TRUE,
                              include.loglik = TRUE,
                              include.deviance = TRUE,
                              include.aic = TRUE,
                              include.bic = TRUE,
                              ...) {
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{negbinirr} objects
#'
#' \code{\link{extract}} method for \code{negbinirr} objects created by the
#' \code{\link[mfx]{negbinirr}} function in the \pkg{mfx} package.
#'
#' @param model A statistical model object.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.deviance Report the deviance?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract negbinirr
#' @aliases extract.negbinirr
#' @export
setMethod("extract", signature = className("negbinirr", "mfx"),
          definition = extract.negbinirr)


# -- extract.poissonirr (mfx) --------------------------------------------------

#' @noRd
extract.poissonirr <- extract.negbinirr


#' \code{\link{extract}} method for \code{poissonirr} objects
#'
#' \code{\link{extract}} method for \code{poissonirr} objects created by the
#' \code{\link[mfx]{poissonirr}} function in the \pkg{mfx} package.
#'
#' @inheritParams extract,negbinirr-method
#'
#' @method extract poissonirr
#' @aliases extract.poissonirr
#' @export
setMethod("extract", signature = className("poissonirr", "mfx"),
          definition = extract.poissonirr)


# -- extract.negbinmfx (mfx) ---------------------------------------------------

#' @noRd
extract.negbinmfx <- function(model,
                              include.nobs = TRUE,
                              include.loglik = TRUE,
                              include.deviance = TRUE,
                              include.aic = TRUE,
                              include.bic = TRUE,
                              ...) {
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{negbinmfx} objects
#'
#' \code{\link{extract}} method for \code{negbinmfx} objects created by the
#' \code{\link[mfx]{negbinmfx}} function in the \pkg{mfx} package.
#'
#' @inheritParams extract,negbinirr-method
#'
#' @method extract negbinmfx
#' @aliases extract.negbinmfx
#' @export
setMethod("extract", signature = className("negbinmfx", "mfx"),
          definition = extract.negbinmfx)


# -- extract.poissonmfx (mfx) --------------------------------------------------

#' @noRd
extract.poissonmfx <- extract.negbinmfx

#' \code{\link{extract}} method for \code{poissonmfx} objects
#'
#' \code{\link{extract}} method for \code{poissonmfx} objects created by the
#' \code{\link[mfx]{poissonmfx}} function in the \pkg{mfx} package.
#'
#' @inheritParams extract,negbinirr-method
#'
#' @method extract poissonmfx
#' @aliases extract.poissonmfx
#' @export
setMethod("extract", signature = className("poissonmfx", "mfx"),
          definition = extract.poissonmfx)


# -- extract.mtergm (btergm) ---------------------------------------------------

#' @noRd
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{mtergm} objects
#'
#' \code{\link{extract}} method for \code{mtergm} objects created by the
#' \code{\link[btergm:btergm]{mtergm}} function in the \pkg{btergm} package.
#'
#' @param model A statistical model object.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract mtergm
#' @aliases extract.mtergm
#' @export
setMethod("extract", signature = className("mtergm", "btergm"),
          definition = extract.mtergm)


# -- extract.nlme (nlme) -------------------------------------------------------

#' @noRd
extract.netlogit <- function(model,
                             include.aic = TRUE,
                             include.bic = TRUE,
                             include.deviance = TRUE,
                             include.nobs = TRUE,
                             ...) {

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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{netlogit} objects
#'
#' \code{\link{extract}} method for \code{netlogit} objects created by the
#' \code{netlogit} function in the \pkg{sna} package.
#'
#' @param model A statistical model object.
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.deviance Report the deviance?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract netlogit
#' @aliases extract.netlogit
#' @export
setMethod("extract", signature = className("netlogit", "sna"),
          definition = extract.netlogit)


# -- extract.nlme (nlme) -------------------------------------------------------

#' @noRd
extract.nlme <- extract.lme

#' \code{\link{extract}} method for \code{nlme} objects
#'
#' \code{\link{extract}} method for \code{nlme} objects created by the
#' \code{\link[nlme]{nlme}} function in the \pkg{nlme} package.
#'
#' @inheritParams extract,lme-method
#'
#' @method extract nlme
#' @aliases extract.nlme
#' @export
setMethod("extract", signature = className("nlme", "nlme"),
          definition = extract.nlme)


# -- extract.oglmx (oglmx) ----------------------------------------------------

#' @noRd
extract.oglmx <- function(model,
                          include.aic = TRUE,
                          include.iterations = TRUE,
                          include.loglik = TRUE,
                          include.nobs = TRUE,
                          include.rsquared = TRUE,
                          ...) {
  s <- summary(model, ...)

  coefficient.names <- names(s$coefficients)
  coefficients <- s$estimate[,1]
  standard.errors <- s$estimate[,2]
  significance <- s$estimate[,4]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    aic <- s$AIC[1]
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    lik <- s$loglikelihood[1]
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  if (include.nobs == TRUE) {
    n <- nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }

  if (include.iterations == TRUE) {
    it <- s$no.iterations
    gof <- c(gof, it)
    gof.names <- c(gof.names, "Iterations")
    gof.decimal <- c(gof.decimal, FALSE)
  }

  if (include.rsquared == TRUE) {
    r2 <- s$McFaddensR2[1]
    gof <- c(gof, r2)
    gof.names <- c(gof.names, "McFadden's R$^2$")
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

#' \code{\link{extract}} method for \code{oglmx} objects
#'
#' \code{\link{extract}} method for \code{oglmx} objects created by the
#' \code{\link[oglmx]{oglmx}} function in the \pkg{oglmx} package.
#'
#' @param model A statistical model object.
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.iterations Report the number of iterations?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.rsquared Report R^2 in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract oglmx
#' @aliases extract.oglmx
#' @export
setMethod("extract", signature = className("oglmx", "oglmx"),
          definition = extract.oglmx)


# -- extract.ols (rms) ---------------------------------------------------------

#' @noRd
extract.ols <- function(model,
                        include.nobs = TRUE,
                        include.rsquared = TRUE,
                        include.adjrs = TRUE,
                        include.fstatistic = FALSE,
                        include.lr = TRUE,
                        ...) {

  names <- attributes(model$coef)$names
  co <- model$coef
  se <- sqrt(diag(model$var))
  pval <- pnorm(abs(model$coef / sqrt(diag(model$var))), lower.tail = FALSE) * 2

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    n <- nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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
    gof.names <- c(gof.names, "Adj. R$^2$")
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

  tr <- createTexreg(coef.names = names,
                     coef = co,
                     se = se,
                     pvalues = pval,
                     gof.names = gof.names,
                     gof = gof,
                     gof.decimal = gof.decimal)
  return(tr)
}

#' \code{\link{extract}} method for \code{ols} objects
#'
#' \code{\link{extract}} method for \code{ols} objects created by the
#' \code{\link[rms]{ols}} function in the \pkg{rms} package.
#'
#' @param model A statistical model object.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.rsquared Report R^2 in the GOF block?
#' @param include.adjrs Report adjusted R^2 in the GOF block?
#' @param include.fstatistic Report the F-statistic in the GOF block?
#' @param include.lr Report likelihood ratio test?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract ols
#' @aliases extract.ols
#' @importFrom stats pnorm summary.lm
#' @export
setMethod("extract", signature = className("ols", "rms"),
          definition = extract.ols)


# -- extract.pcce (plm) --------------------------------------------------------

#' @noRd
extract.pcce <- function(model,
                         include.r.squared = TRUE,
                         include.sumsquares = TRUE,
                         include.nobs = TRUE,
                         ...) {

  s <- summary(model, ...)
  coefnames <- rownames(s$CoefTable)
  co <- s$CoefTable[, 1]
  se <- s$CoefTable[, 2]
  pval <- s$CoefTable[, 4]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (isTRUE(include.r.squared)) {
    rs <- s$r.squared
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "HPY R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (isTRUE(include.sumsquares)) {
    gof <- c(gof, s$ssr, s$tss)
    gof.names <- c(gof.names, "RSS", "TSS")
    gof.decimal <- c(gof.decimal, TRUE, TRUE)
  }
  if (isTRUE(include.nobs)) {
    n <- nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }

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
}

#' \code{\link{extract}} method for \code{pcce} objects
#'
#' \code{\link{extract}} method for \code{pcce} objects created by the
#' \code{\link[plm]{pcce}} function in the \pkg{plm} package.
#'
#' @param model A statistical model object.
#' @param include.r.squared Report the HPY R-squared statistic in the GOF block?
#' @param include.sumsquares Report the total and residual sum of squares in the
#'   GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract pcce
#' @aliases extract.pcce
#' @export
setMethod("extract", signature = className("pcce", "plm"),
          definition = extract.pcce)


# -- extract.pglm (pglm) -------------------------------------------------------

#' @noRd
extract.pglm <- function(model,
                         include.aic = TRUE,
                         include.loglik = TRUE,
                         include.nobs = TRUE,
                         ...) {
  s <- summary(model, ...)

  coefficient.names <- rownames(s$estimate)
  coefficients <- s$estimate[, 1]
  standard.errors <- s$estimate[, 2]
  significance <- s$estimate[, 4]

  aic <- AIC(model)
  lik <- logLik(model)[1]
  n <- length(model$gradientObs[,1])

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{pglm} objects
#'
#' \code{\link{extract}} method for \code{pglm} objects created by the
#' \code{\link[pglm]{pglm}} function in the \pkg{pglm} package.
#'
#' @param model A statistical model object.
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract pglm
#' @aliases extract.pglm
#' @export
setMethod("extract", signature = className("pglm", "pglm"),
          definition = extract.pglm)


# -- extract.pgmm (plm) --------------------------------------------------------

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
    gof.names <- c(gof.names, "n", "T", "Num. obs.", "Num. obs. used")
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

#' \code{\link{extract}} method for \code{pgmm} objects
#'
#' \code{\link{extract}} method for \code{pgmm} objects created by the
#' \code{\link[plm]{pgmm}} function in the \pkg{plm} package.
#'
#' @param model A statistical model object.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.sargan Report the Sargan test?
#' @param include.wald Report the Wald statistic?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract pgmm
#' @aliases extract.pgmm
#' @export
setMethod("extract", signature = className("pgmm", "plm"),
          definition = extract.pgmm)


# -- extract.plm (plm) ---------------------------------------------------------

#' @noRd
extract.plm <- function(model,
                        include.rsquared = TRUE,
                        include.adjrs = TRUE,
                        include.nobs = TRUE,
                        include.variance = TRUE,
                        ...) {
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
    gof.names <- c(gof.names, "Adj. R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{plm} objects
#'
#' \code{\link{extract}} method for \code{plm} objects created by the
#' \code{\link[plm]{plm}} function in the \pkg{plm} package.
#'
#' @param model A statistical model object.
#' @param include.rsquared Report R^2 in the GOF block?
#' @param include.adjrs Report adjusted R^2 in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.variance Report group variances?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract plm
#' @aliases extract.plm
#' @export
setMethod("extract", signature = className("plm", "plm"),
          definition = extract.plm)


# -- extract.pmg (plm) ----------------------------------------------------

#' @noRd
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{pmg} objects
#'
#' \code{\link{extract}} method for \code{pmg} objects created by the
#' \code{\link[plm]{pmg}} function in the \pkg{plm} package.
#'
#' @param model A statistical model object.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract pmg
#' @aliases extract.pmg
#' @importFrom stats pt
#' @export
setMethod("extract", signature = className("pmg", "plm"),
          definition = extract.pmg)


# -- extract.polr (MASS) -------------------------------------------------------

#' @noRd
extract.polr <- function(model,
                         include.thresholds = TRUE,
                         include.aic = TRUE,
                         include.bic = TRUE,
                         include.loglik = TRUE,
                         include.deviance = TRUE,
                         include.nobs = TRUE,
                         ...) {
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{polr} objects
#'
#' \code{\link{extract}} method for \code{polr} objects created by the
#' \code{\link[MASS]{polr}} function in the \pkg{MASS} package.
#'
#' @param model A statistical model object.
#' @param include.thresholds Report thresholds in the GOF block?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.deviance Report the deviance?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract polr
#' @aliases extract.polr
#' @importFrom stats pnorm AIC BIC logLik deviance
#' @export
setMethod("extract", signature = className("polr", "MASS"),
          definition = extract.polr)


# -- extract.prais (prais) -----------------------------------------------------

#' @noRd
extract.prais <- function(model,
                          include.rsquared = TRUE,
                          include.adjrs = TRUE,
                          include.nobs = TRUE,
                          ...) {
  s <- summary(model, ...)

  rho <- mean(s$rho[nrow(s$rho), ])

  coefficient.names <- c(rownames(coef(s)), "rho")
  coefficients <- c(coef(s)[, 1], rho)
  standard.errors <- c(coef(s)[, 2], NA_real_)
  significance <- c(coef(s)[, 4], NA_real_)

  rs <- s$r.squared
  adj <- s$adj.r.squared
  n <- length(model$residuals)

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
    gof.names <- c(gof.names, "Adj. R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{prais} objects
#'
#' \code{\link{extract}} method for \code{prais} objects created by the
#' \code{\link[prais]{prais_winsten}} function in the \pkg{prais} package.
#'
#' @param model A statistical model object.
#' @param include.rsquared Report R^2 in the GOF block?
#' @param include.adjrs Report adjusted R^2 in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract prais
#' @aliases extract.prais
#' @export
setMethod("extract", signature = className("prais", "prais"),
          definition = extract.prais)


# -- extract.rem.dyad (relevent) -----------------------------------------------

# extension for rem.dyad objects (relevent package)
extract.rem.dyad <- function(model,
                             include.nvertices = TRUE,
                             include.events = TRUE,
                             include.aic = TRUE,
                             include.aicc = TRUE,
                             include.bic = TRUE,
                             ...) {

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
    gof.names <- c(gof.names, "Num. nodes")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.events == TRUE) {
    num.events <- model$m
    gof <- c(gof, num.events)
    gof.names <- c(gof.names, "Num. events")
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

#' \code{\link{extract}} method for \code{rem.dyad} objects
#'
#' \code{\link{extract}} method for \code{rem.dyad} objects created by the
#' \code{\link[relevent]{rem.dyad}} function in the \pkg{relevent} package.
#'
#' @param model A statistical model object.
#' @param include.nvertices Report the number of vertices in a STERGM?
#' @param include.events Report the number of events in the GOF block?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.aicc Report AICC in the GOF block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract rem.dyad
#' @aliases extract.rem.dyad
#' @importFrom stats pnorm
#' @export
setMethod("extract", signature = className("rem.dyad", "relevent"),
          definition = extract.rem.dyad)


# -- extract.remstimate (remstimate) -------------------------------------------

#' @noRd
extract.remstimate <- function(model, include.nobs = TRUE, include.aic = TRUE,
                               include.aicc = TRUE, include.bic = TRUE,
                               include.loglik = TRUE, post.mean = FALSE, ...) {
  op <- capture.output(s <- summary(model))
  if (length(model) > 2) {
    models <- list(model)
  } else {
    models <- model
  }
  tr <- list()
  for (i in 1:length(models)) {
    coefs <- models[[i]]$coefficients
    coefnames <- attributes(models[[i]]$coefficients)$names

    if (length(model) > 2) {
      coefstab <- s$coefsTab
    } else {
      coefstab <- s$coefsTab[[names(models)[[i]]]]
    }
    se <- coefstab[, 2]

    gof <- numeric()
    gof.names <- character()
    gof.decimal <- logical()
    if (include.nobs == TRUE) {
      gof <- c(gof, models[[i]]$df.null)
      gof.names <- c(gof.names, "Num. obs.")
      gof.decimal <- c(gof.decimal, FALSE)
    }
    if (attributes(model)$method[1] %in% c("MLE", "GDADAMAX")) {
      pval <- coefstab[, 4]
      if (include.aic == TRUE) {
        gof <- c(gof, models[[i]]$AIC)
        gof.names <- c(gof.names, "AIC")
        gof.decimal <- c(gof.decimal, TRUE)
      }
      if (include.aicc == TRUE) {
        gof <- c(gof, models[[i]]$AICC)
        gof.names <- c(gof.names, "AICc")
        gof.decimal <- c(gof.decimal, TRUE)
      }
      if (include.bic == TRUE) {
        gof <- c(gof, models[[i]]$BIC)
        gof.names <- c(gof.names, "BIC")
        gof.decimal <- c(gof.decimal, TRUE)
      }
    } else if (attributes(model)$method[1] %in% c("HMC", "BSIR")) {
      pval <- coefstab[, 6]
      ci.l <- coefstab[, 3]
      ci.u <- coefstab[, 5]
      if (post.mean == TRUE) {
        coefs <- coefstab[, 4]
      }
    }
    if (include.loglik == TRUE) {
      gof <- c(gof, models[[i]]$loglik)
      gof.names <- c(gof.names, "Log Likelihood")
      gof.decimal <- c(gof.decimal, TRUE)
    }

    if (length(model) > 2) {
      tr[[i]] <- createTexreg(
        coef.names = coefnames,
        coef = coefs,
        se = se,
        pvalues = pval,
        gof.names = gof.names,
        gof = gof,
        gof.decimal = gof.decimal
      )
    } else {
      tr[[i]] <- createTexreg(
        coef.names = coefnames,
        coef = coefs,
        se = se,
        pvalues = pval,
        gof.names = gof.names,
        gof = gof,
        gof.decimal = gof.decimal,
        model.name = names(model)[i]
      )
    }
  }
  if (length(tr) == 1) {
    return(tr[[1]])
  } else {
    return(tr)
  }
}

#' \code{\link{extract}} method for \code{remstimate} objects
#'
#' \code{\link{extract}} method for \code{remstimate} objects created by the
#' \code{\link[remstimate]{remstimate}} function in the \pkg{remstimate}
#' package.
#'
#' @param model A statistical model object.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.aicc Report Corrected Akaike's Information Criterion (AICc) in
#'   the GOF block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param post.mean Report the posterior means, rather than the posterior modes,
#'   as coefficients?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract remstimate
#' @aliases extract.remstimate
#' @importFrom utils capture.output
#' @export
setMethod("extract", signature = className("remstimate", "remstimate"),
          definition = extract.remstimate)


# -- extract.rlm (MASS) --------------------------------------------------------

#' @noRd
extract.rlm <- function (model, include.nobs = TRUE, ...) {
  s <- summary(model, ...)

  names <- rownames(s$coef)
  co <- s$coef[, 1]
  se <- s$coef[, 2]
  pval <- 1 - pchisq((co / se)^2, 1)

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    n <- nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{rlm} objects
#'
#' \code{\link{extract}} method for \code{rlm} objects created by the
#' \code{\link[MASS]{rlm}} function in the \pkg{MASS} package.
#'
#' @param model A statistical model object.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract rlm
#' @aliases extract.rlm
#' @importFrom stats pchisq
#' @export
setMethod("extract", signature = className("rlm", "MASS"),
          definition = extract.rlm)


# -- extract.rq (quantreg) -----------------------------------------------------

#' @noRd
extract.rq <- function(model,
                       include.nobs = TRUE,
                       include.percentile = TRUE,
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{rq} objects
#'
#' \code{\link{extract}} method for \code{rq} objects created by the
#' \code{rq} function in the \pkg{quantreg} package.
#'
#' @param model A statistical model object.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.percentile Report the percentile (tau)?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract rq
#' @aliases extract.rq
#' @export
setMethod("extract", signature = className("rq", "quantreg"),
          definition = extract.rq)


# -- extract.Sarlm (spatialreg) ------------------------------------------------

#' @noRd
extract.Sarlm <- function(model,
                          include.nobs = TRUE,
                          include.loglik = TRUE,
                          include.aic = TRUE,
                          include.lr = TRUE,
                          include.wald = TRUE,
                          ...) {
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
    gof.names <- c(gof.names, "Num. obs.", "Parameters")
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

#' \code{\link{extract}} method for \code{Sarlm} objects
#'
#' \code{\link{extract}} method for \code{Sarlm} objects created by the
#' \code{\link[spatialreg:ML_models]{lagsarlm}} function in the \pkg{spatialreg}
#' package.
#'
#' @param model A statistical model object.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.lr Report likelihood ratio test?
#' @param include.wald Report the Wald statistic?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract Sarlm
#' @aliases extract.Sarlm
#' @export
#' @importFrom stats pnorm AIC
setMethod("extract", signature = className("Sarlm", "spatialreg"),
          definition = extract.Sarlm)


# -- extract.sclm (ordinal) ----------------------------------------------------

#' @noRd
extract.sclm <- extract.clm

#' \code{\link{extract}} method for \code{sclm} objects
#'
#' \code{\link{extract}} method for \code{sclm} objects created by the
#' \code{\link[ordinal]{clm}} function in the \pkg{ordinal} package.
#'
#' @inheritParams extract,clm-method
#'
#' @method extract sclm
#' @aliases extract.sclm
#' @export
setMethod("extract", signature = className("sclm", "ordinal"),
          definition = extract.sclm)


# -- extract.selection (sampleSelection) ---------------------------------------

#' @noRd
extract.selection <- function(model,
                              prefix = TRUE,
                              include.selection = TRUE,
                              include.outcome = TRUE,
                              include.errors = TRUE,
                              include.aic = TRUE,
                              include.bic = TRUE,
                              include.loglik = TRUE,
                              include.rsquared = TRUE,
                              include.adjrs = TRUE,
                              include.nobs = TRUE,
                              ...) {

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
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.rsquared == TRUE && "rSquared" %in% names(s)) {
    gof <- c(gof, s$rSquared$R2)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.adjrs == TRUE && "rSquared" %in% names(s)) {
    gof <- c(gof, s$rSquared$R2adj)
    gof.names <- c(gof.names, "Adj. R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    if(model$tobitType == 5) {
      gof <- c(gof, s$param$nObs, s$param$N1, s$param$N2)
    } else if(model$tobitType == 2) {
      gof <- c(gof, s$param$nObs, s$param$N0, s$param$N1)
    }
    gof.names <- c(gof.names, "Num. obs.", "Censored", "Observed")
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

#' \code{\link{extract}} method for \code{selection} objects
#'
#' \code{\link{extract}} method for \code{selection} objects created by the
#' \code{selection} function in the \pkg{sampleSelection} package.
#'
#' @param model A statistical model object.
#' @param prefix Include prefix before the label of the coefficient in order to
#'   identify the current model component?
#' @param include.selection Report the selection component of a sample selection
#'   model?
#' @param include.outcome Report the outcome component of a sample selection
#'   model?
#' @param include.errors Report the error terms of a sample selection model?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.rsquared Report R^2 in the GOF block?
#' @param include.adjrs Report adjusted R^2 in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#'   block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract selection
#' @aliases extract.selection
#' @export
setMethod("extract", signature = className("selection", "sampleSelection"),
          definition = extract.selection)


# -- extract.sienaFit (RSiena) -------------------------------------------------

#' @noRd
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

#' \code{\link{extract}} method for \code{sienaFit} objects
#'
#' \code{\link{extract}} method for \code{sienaFit} objects created by the
#' \code{siena07} function in the \pkg{RSiena} package.
#'
#' @param model A statistical model object.
#' @param include.iterations Report the number of iterations?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract sienaFit
#' @aliases extract.sienaFit
#' @importFrom stats pnorm
#' @export
setMethod("extract", signature = className("sienaFit", "RSiena"),
          definition = extract.sienaFit)


# -- extract.simex (simex) -----------------------------------------------------

#' @noRd
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{simex} objects
#'
#' \code{\link{extract}} method for \code{simex} objects created by the
#' \code{\link[simex]{simex}} function in the \pkg{simex} package.
#'
#' @param model A statistical model object.
#' @param jackknife Use Jackknife variance instead of asymptotic variance?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract simex
#' @aliases extract.simex
#' @export
setMethod("extract", signature = className("simex", "simex"),
          definition = extract.simex)


# -- extract.speedglm (speedglm) -----------------------------------------------

#' @noRd
extract.speedglm <- function(model, include.aic = TRUE, include.bic = TRUE,
                             include.loglik = TRUE, include.deviance = TRUE,
                             include.nobs = TRUE, ...) {
  s <- summary(model, ...)
  coefficient.names <- rownames(s$coefficients)
  coefficients <- model$coefficients
  standard.errors <- s$coefficients[, 2]
  tval <- model$coefficients / standard.errors
  if (model$family$family %in% c("binomial", "poisson")) {
    significance <- 2 * pnorm(abs(tval), lower.tail = FALSE)
  } else {
    significance <- 2 * pt(abs(tval), df = model$df, lower.tail = FALSE)
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
    dev <- deviance(model)
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    n <- nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{speedglm} objects
#'
#' \code{\link{extract}} method for \code{speedglm} objects created by the
#' \code{\link[speedglm]{speedglm}} function in the \pkg{speedglm}
#' package.
#'
#' @inheritParams extract,glm-method
#'
#' @method extract speedglm
#' @aliases extract.speedglm
#' @export
setMethod("extract",  signature = className("speedglm", "speedglm"),
          definition = extract.speedglm)


# -- extract.speedlm (speedglm) -----------------------------------------------

#' @noRd
extract.speedlm <- extract.lm

#' \code{\link{extract}} method for \code{speedlm} objects
#'
#' \code{\link{extract}} method for \code{speedlm} objects created by the
#' \code{\link[speedglm]{speedlm}} function in the \pkg{speedglm}
#' package.
#'
#' @inheritParams extract,lm-method
#'
#' @method extract speedlm
#' @aliases extract.speedlm
#' @export
setMethod("extract",  signature = className("speedlm", "speedglm"),
          definition = extract.speedlm)


# -- extract.stergm (tergm) ----------------------------------------------------

#' @noRd
extract.stergm <- function(model,
                           beside = FALSE,
                           include.formation = TRUE,
                           include.dissolution = TRUE,
                           include.nvertices = TRUE,
                           include.aic = FALSE,
                           include.bic = FALSE,
                           include.loglik = FALSE,
                           ...) {
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
      gof.names <- c(gof.names, "Num. vertices")
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
      f.gof.names <- c(f.gof.names, "Num. vertices")
      f.gof.decimal <- c(f.gof.decimal, FALSE)
      d.gof <- c(d.gof, nvertices)
      d.gof.names <- c(d.gof.names, "Num. vertices")
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

#' \code{\link{extract}} method for \code{stergm} objects
#'
#' \code{\link{extract}} method for \code{stergm} objects created by the
#' \code{stergm} function in the \pkg{tergm} package.
#'
#' @param model A statistical model object.
#' @param beside Arrange the model terms below each other or beside each other?
#'   In a \code{stergm} model, the formation and dissolution coefficients can be
#'   arranged in two columns of the table.
#' @param include.formation Report the coefficients for the formation process in
#'   a STERGM?
#' @param include.dissolution Report the coefficients for the dissolution
#'   process in a STERGM?
#' @param include.nvertices Report the number of vertices in a STERGM?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract stergm
#' @aliases extract.stergm
#' @export
setMethod("extract", signature = className("stergm", "tergm"),
          definition = extract.stergm)


# -- extract.summary.lm (stats) ------------------------------------------------

#' @noRd
extract.summary.lm <- function (model,
                                include.rsquared = TRUE,
                                include.adjrs = TRUE,
                                include.nobs = TRUE,
                                include.fstatistic = FALSE,
                                include.rmse = TRUE,
                                ...) {
  s <- model
  names <- rownames(s$coef)
  co <- s$coef[, 1]
  se <- s$coef[, 2]
  pval <- s$coef[, 4]

  rs <- s$r.squared
  adj <- s$adj.r.squared
  n <- length(s$residuals)

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
    gof.names <- c(gof.names, "Adj. R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{summary.lm} objects
#'
#' \code{\link{extract}} method for \code{summary.lm} objects created by the
#' \code{summary} method for \code{lm} objects, defined in the \pkg{stats}
#' package (see \code{\link[stats]{summary.lm}}).
#'
#' @param model A statistical model object.
#' @param include.rsquared Report R^2 in the GOF block?
#' @param include.adjrs Report adjusted R^2 in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.fstatistic Report the F-statistic in the GOF block?
#' @param include.rmse Report the root mean square error (RMSE; = residual
#'   standard deviation) in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract summary.lm
#' @aliases extract.summary.lm
#' @export
setMethod("extract",  signature = className("summary.lm", "stats"),
          definition = extract.summary.lm)


# -- extract.survreg (survival) ------------------------------------------------

#' @noRd
extract.survreg <- function(model,
                            include.aic = TRUE,
                            include.bic = TRUE,
                            include.loglik = TRUE,
                            include.deviance = TRUE,
                            include.nobs = TRUE,
                            ...) {

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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{survreg} objects
#'
#' \code{\link{extract}} method for \code{survreg} objects created by the
#' \code{\link[survival]{survreg}} function in the \pkg{survival} package.
#'
#' @param model A statistical model object.
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.deviance Report the deviance?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract survreg
#' @aliases extract.survreg
#' @importFrom stats AIC BIC
#' @export
setMethod("extract", signature = className("survreg", "survival"),
          definition = extract.survreg)


# -- extract.survreg.penal (survival) ------------------------------------------

#' @noRd
extract.survreg.penal <- extract.survreg

#' \code{\link{extract}} method for \code{survreg.penal} objects
#'
#' \code{\link{extract}} method for \code{survreg.penal} objects created by the
#' \code{\link[survival]{survreg}} function in the \pkg{survival} package.
#'
#' @inheritParams extract,survreg-method
#'
#' @method extract survreg.penal
#' @aliases extract.survreg.penal
#' @export
setMethod("extract", signature = className("survreg.penal", "survival"),
          definition = extract.survreg.penal)


# -- extract.svyglm (survey) ---------------------------------------------------

#' @noRd
extract.svyglm <- function(model,
                           include.aic = FALSE,
                           include.bic = FALSE,
                           include.loglik = FALSE,
                           include.deviance = TRUE,
                           include.dispersion = TRUE,
                           include.nobs = TRUE,
                           ...) {
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{svyglm} objects
#'
#' \code{\link{extract}} method for \code{svyglm} objects created by the
#' \code{svyglm} function in the \pkg{survey} package.
#'
#' @param model A statistical model object.
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.deviance Report the deviance?
#' @param include.dispersion Report the dispersion parameter?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract svyglm
#' @aliases extract.svyglm
#' @export
setMethod("extract", signature = className("svyglm", "survey"),
          definition = extract.svyglm)


# -- extract.systemfit (systemfit) ---------------------------------------------------

#' @noRd
extract.systemfit <- function(model,
                              include.rsquared = TRUE,
                              include.adjrs = TRUE,
                              include.nobs = TRUE,
                              beside = FALSE,
                              include.suffix = FALSE,
                              ...) {
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
        gof.names <- c(gof.names, "Adj. R$^2$")
        gof.decimal <- c(gof.decimal, TRUE)
      }
      if (include.nobs == TRUE) {
        n <- length(s$residuals)  # extract number of observations
        gof <- c(gof, n)
        gof.names <- c(gof.names, "Num. obs.")
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
            paste0("Adj. R$^2$ (", eq$eqnLabel, ")"),
            paste0(eq$eqnLabel, ": Adj. R$^2$")))
        gof.decimal <- c(gof.decimal, TRUE)
      }
    }

    if (include.nobs == TRUE) {  # number of observations
      gof <- c(gof, nobs(model))
      gof.names <- c(gof.names, "Num. obs. (total)")
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

#' \code{\link{extract}} method for \code{systemfit} objects
#'
#' \code{\link{extract}} method for \code{systemfit} objects created by the
#' \code{systemfit} function in the \pkg{systemfit} package.
#'
#' @param model A statistical model object.
#' @param include.rsquared Report R^2 in the GOF block?
#' @param include.adjrs Report adjusted R^2 in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param beside Arrange the model terms below each other or beside each other,
#'   in separate columns?
#' @param include.suffix Report the name of the current model in parentheses
#'   after each model term (instead of before the model term)?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract systemfit
#' @aliases extract.systemfit
#' @export
setMethod("extract", signature = className("systemfit", "systemfit"),
          definition = extract.systemfit)


# -- extract.texreg (texreg) ---------------------------------------------------

#' @noRd
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
    gof.decimal = model@gof.decimal,
    model.name = model@model.name
  )
  return(tr)
}

#' \code{\link{extract}} method for \code{texreg} objects
#'
#' \code{\link{extract}} method for \code{texreg} objects created by the
#' \code{\link[texreg]{extract}} function in the \pkg{texreg} package.
#'
#' @param model A statistical model object.
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract texreg
#' @aliases extract.texreg
#' @export
setMethod("extract", signature = className("texreg", "texreg"),
          definition = extract.texreg)


# -- extract.tobit (AER) -------------------------------------------------------

#' @noRd
extract.tobit <- function(model,
                          include.aic = TRUE,
                          include.bic = TRUE,
                          include.loglik = TRUE,
                          include.deviance = TRUE,
                          include.nobs = FALSE,
                          include.censnobs = TRUE,
                          include.wald = TRUE,
                          ...) {
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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{tobit} objects
#'
#' \code{\link{extract}} method for \code{tobit} objects created by the
#' \code{\link[AER]{tobit}} function in the \pkg{AER} package.
#'
#' @param model A statistical model object.
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.deviance Report the deviance?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.censnobs Report the total, right-censored, left-censored, and
#'   uncensored number of observations?
#' @param include.wald Report the Wald statistic?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract tobit
#' @aliases extract.tobit
#' @importFrom stats AIC BIC logLik deviance nobs
#' @export
setMethod("extract", signature = className("tobit", "AER"),
          definition = extract.tobit)


# -- extract.truncreg (truncreg) -----------------------------------------------

#' @noRd
extract.truncreg <- function(model,
                             include.nobs = TRUE,
                             include.loglik = TRUE,
                             include.aic = TRUE,
                             include.bic = TRUE,
                             ...) {
  s <- summary(model, ...)
  coef_names <- rownames(s$coefficients)
  coefs <- s$coefficients[, 1]
  se <- s$coefficients[, 2]
  pval <- s$coefficients[, 4]

  gof_names <- character()
  gof_stats <- numeric()
  gof_decimal <- logical()
  if (isTRUE(include.nobs)) {
    gof_names <- c(gof_names, "Num. obs.")
    gof_stats <- c(gof_stats, s$nobs)
    gof_decimal <- c(gof_decimal, FALSE)
  }
  if (isTRUE(include.loglik)) {
    gof_names <- c(gof_names, "Log Likelihood")
    gof_stats <- c(gof_stats, s$logLik)
    gof_decimal <- c(gof_decimal, TRUE)
  }
  if (isTRUE(include.aic)) {
    gof_names <- c(gof_names, "AIC")
    aic <- (-2 * s$logLik) + (2 * length(model$coefficients))
    gof_stats <- c(gof_stats, aic)
    gof_decimal <- c(gof_decimal, TRUE)
  }
  if (isTRUE(include.bic)) {
    gof_names <- c(gof_names, "BIC")
    bic <- (-2 * s$logLik) + (log(model$nobs) * length(model$coefficients))
    gof_stats <- c(gof_stats, bic)
    gof_decimal <- c(gof_decimal, TRUE)
  }

  createTexreg(coef.names = coef_names,
               coef = coefs,
               se = se,
               pvalues = pval,
               gof.names = gof_names,
               gof = gof_stats,
               gof.decimal = gof_decimal)
}

#' \code{\link{extract}} method for \code{truncreg} objects
#'
#' \code{\link{extract}} method for \code{truncreg} objects created by the
#' \code{\link[truncreg]{truncreg}} function in the \pkg{truncreg} package.
#'
#' @param model A statistical model object.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract truncreg
#' @aliases extract.truncreg
#' @export
setMethod("extract", signature = className("truncreg", "truncreg"),
          definition = extract.truncreg)


# -- extract.vglm (VGAM) -------------------------------------------------------

#' @noRd
extract.vglm <- function(model,
                         include.loglik = TRUE,
                         include.df = TRUE,
                         include.nobs = TRUE,
                         ...) {

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
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{vglm} objects
#'
#' \code{\link{extract}} method for \code{vglm} objects created by the
#' \code{\link[VGAM]{vglm}} function in the \pkg{VGAM} package.
#'
#' @param model A statistical model object.
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.df Report the degrees of freedom?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract vglm
#' @aliases extract.vglm
#' @author Christoph Riedl <c.riedl@neu.edu>
#' @export
setMethod("extract", signature = className("vglm", "VGAM"),
          definition = extract.vglm)


# -- extract.weibreg (eha) -----------------------------------------------------

#' @noRd
extract.weibreg <- function(model,
                            include.aic = TRUE,
                            include.loglik = TRUE,
                            include.lr = TRUE,
                            include.nobs = TRUE,
                            include.events = TRUE,
                            include.trisk = TRUE,
                            ...) {

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
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    n <- nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.events == TRUE) {
    ev <- model$n.events
    gof <- c(gof, ev)
    gof.names <- c(gof.names, "Num. events")
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

#' \code{\link{extract}} method for \code{weibreg} objects
#'
#' \code{\link{extract}} method for \code{weibreg} objects created by the
#' \code{\link[eha]{weibreg}} function in the \pkg{eha} package.
#'
#' @param model A statistical model object.
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.lr Report likelihood ratio test?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.events Report the number of events in the GOF block?
#' @param include.trisk Report the total time at risk (in event-history models)?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract weibreg
#' @aliases extract.weibreg
#' @importFrom stats pchisq nobs
#' @export
setMethod("extract", signature = className("weibreg", "eha"),
          definition = extract.weibreg)


# -- extract.phreg (eha) -------------------------------------------------------

#' @noRd
extract.phreg <- extract.weibreg

#' \code{\link{extract}} method for \code{phreg} objects
#'
#' \code{\link{extract}} method for \code{phreg} objects created by the
#' \code{\link[eha]{phreg}} function in the \pkg{eha} package.
#'
#' @inheritParams extract,weibreg-method
#'
#' @method extract phreg
#' @aliases extract.phreg
#' @importFrom stats pchisq nobs
#' @export
setMethod("extract", signature = className("phreg", "eha"),
          definition = extract.phreg)


# -- extract.aftreg (eha) ------------------------------------------------------

#' @noRd
extract.aftreg <- extract.weibreg

#' \code{\link{extract}} method for \code{aftreg} objects
#'
#' \code{\link{extract}} method for \code{aftreg} objects created by the
#' \code{\link[eha]{aftreg}} function in the \pkg{eha} package.
#'
#' @inheritParams extract,weibreg-method
#'
#' @method extract aftreg
#' @aliases extract.aftreg
#' @importFrom stats pchisq nobs
#' @export
setMethod("extract", signature = className("aftreg", "eha"),
          definition = extract.aftreg)


# -- extract.coxreg (eha) ------------------------------------------------------

#' @noRd
extract.coxreg <- extract.weibreg

#' \code{\link{extract}} method for \code{coxreg} objects
#'
#' \code{\link{extract}} method for \code{coxreg} objects created by the
#' \code{\link[eha]{coxreg}} function in the \pkg{eha} package.
#'
#' @inheritParams extract,weibreg-method
#'
#' @method extract coxreg
#' @aliases extract.coxreg
#' @importFrom stats pchisq nobs
#' @export
setMethod("extract", signature = className("coxreg", "eha"),
          definition = extract.coxreg)


# -- extract.wls (metaSEM) -----------------------------------------------------

#' @noRd
extract.wls <- function(model,
                        include.statistics = TRUE,
                        include.nobs = TRUE,
                        include.aic = TRUE,
                        include.bic = TRUE,
                        ...) {

  coefnames <- rownames(summary(model)$coef)
  coefs <- summary(model)$coef$Estimate
  se <- as.numeric(summary(model)$coef$`Std.Error`)
  pval <- summary(model)$coef$`Pr(>|z|)`
  ci.l <- summary(model)$coef$lbound
  ci.u <- summary(model)$coef$ubound
  if (all(is.na(se))) {
    se <- numeric(0)
  }
  if (all(is.na(pval))) {
    pval <- numeric(0)
  }
  if (all(is.na(ci.l)) && all(is.na(ci.u))) {
    ci.l <- numeric(0)
    ci.u <- numeric(0)
  }
  # prefer SE and p-values if available
  if (length(se) > 0 && length(se) == length(pval) &&
      (length(ci.u) > 0 || length(ci.l) > 0)) {
    ci.l <- numeric(0)
    ci.u <- numeric(0)
  }

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.statistics == TRUE) {
    chi <- summary(model)$stat["Chi-square of target model", 1]
    if (!is.na(chi) && !is.null(chi) && is.finite(chi)) {
      gof <- c(gof, chi)
      gof.names <- c(gof.names, "Chi-square of target model")
      gof.decimal <- c(gof.decimal, TRUE)
    }
    dfs <- summary(model)$stat["DF of target model", 1]
    if (!is.na(dfs) && !is.null(dfs) && is.finite(dfs)) {
      gof <- c(gof, dfs)
      gof.names <- c(gof.names, "DF of target model")
      gof.decimal <- c(gof.decimal, FALSE)
    }
    chi.pval <- summary(model)$stat["p value of target model", 1]
    if (!is.na(chi.pval) && !is.null(chi.pval) && is.finite(chi.pval)) {
      gof <- c(gof, chi.pval)
      gof.names <- c(gof.names, "Chi-square p-value")
      gof.decimal <- c(gof.decimal, TRUE)
    }
    rmsea <- summary(model)$stat["RMSEA", 1]
    if (!is.na(rmsea) && !is.null(rmsea) && is.finite(rmsea)) {
      gof <- c(gof, rmsea)
      gof.names <- c(gof.names, "RMSEA")
      gof.decimal <- c(gof.decimal, TRUE)
    }
    rmseall  <- summary(model)$stat["RMSEA lower 95% CI", 1]
    if (!is.na(rmseall) && !is.null(rmseall) && is.finite(rmseall)) {
      gof <- c(gof, rmseall)
      gof.names <- c(gof.names, "RMSEA lower 95 percent CI")
      gof.decimal <- c(gof.decimal, TRUE)
    }
    rmseaul  <- summary(model)$stat["RMSEA upper 95% CI", 1]
    if (!is.na(rmseaul) && !is.null(rmseaul) && is.finite(rmseaul)) {
      gof <- c(gof, rmseaul)
      gof.names <- c(gof.names, "RMSEA upper 95 percent CI")
      gof.decimal <- c(gof.decimal, TRUE)
    }
    cfi <- summary(model)$stat["CFI", 1]
    if (!is.na(cfi) && !is.null(cfi) && is.finite(cfi)) {
      gof <- c(gof, cfi)
      gof.names <- c(gof.names, "CFI")
      gof.decimal <- c(gof.decimal, TRUE)
    }
    # Compute average variance extracted
    # Based on: https://cran.r-project.org/web/packages/cSEM/vignettes/Using-assess.html#ave
    # Given that the model uses standardized loadings: Var() = 1 and the denominator reduces to K
    # "Empirically, it is the mean square of the standardized loading"
    mat <- model$mx.fit$impliedS1$result
    if (is.null(mat)) {
      ave <- NULL
    } else {
      ave <- mean(mat[nrow(mat), -ncol(mat)]^2)
    }
    if (!is.null(ave) && !is.na(ave) && is.finite(ave)) {
      gof <- c(gof, ave)
      gof.names <- c(gof.names, "Average variance extracted")
      gof.decimal <- c(gof.decimal, TRUE)
    }
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, summary(model)$stat["Sample size", 1])
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.aic == TRUE) {
    aic <- summary(model)$stat["AIC", 1]
    if (!is.null(aic) && !is.na(aic) && length(aic) == 1 && is.finite(aic)) {
      gof <- c(gof, aic)
      gof.names <- c(gof.names, "AIC")
      gof.decimal <- c(gof.decimal, TRUE)
    }
  }
  if (include.bic == TRUE) {
    bic <- summary(model)$stat["BIC", 1]
    if (!is.null(bic) && !is.na(bic) && length(bic) == 1 && is.finite(bic)) {
      gof <- c(gof, bic)
      gof.names <- c(gof.names, "BIC")
      gof.decimal <- c(gof.decimal, TRUE)
    }
  }

  tr <- createTexreg(
    coef.names = coefnames,
    coef = coefs,
    se = se,
    pvalues = pval,
    ci.low = ci.l,
    ci.up = ci.u,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
  return(tr)
}

#' \code{\link{extract}} method for \code{wls} objects
#'
#' \code{\link{extract}} method for \code{wls} objects created by the
#' \code{\link[metaSEM]{wls}} function in the \pkg{metaSEM} package.
#'
#' @param model A statistical model object.
#' @param include.statistics Report RMSEA and other GOF statistics?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.bic Report BIC?
#' @param include.aic Report AIC?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract wls
#' @aliases extract.wls
#' @author Christoph Riedl <c.riedl@neu.edu>
#' @author Philip Leifeld
#' @export
setMethod("extract", signature = className("wls", "metaSEM"),
          definition = extract.wls)


# -- extract.zeroinfl (pscl) ---------------------------------------------------

#' @noRd
extract.zeroinfl <- function(model,
                             beside = FALSE,
                             include.count = TRUE,
                             include.zero = TRUE,
                             include.aic = TRUE,
                             include.loglik = TRUE,
                             include.nobs = TRUE,
                             ...) {

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
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    n <- s$n
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
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

#' \code{\link{extract}} method for \code{zeroinfl} objects
#'
#' \code{\link{extract}} method for \code{zeroinfl} objects created by the
#' \code{zeroinfl} function in the \pkg{pscl} package.
#'
#' @param model A statistical model object.
#' @param beside Arrange the model terms below each other or beside each other?
#'   The binary model parameters and the count parameters can be displayed in
#'   two separate columns of the table.
#' @param include.count Report the count parameters in the coefficients block
#'   (before the binary part for the zeros)?
#' @param include.zero Should the binary part of the model be included in the
#'   coefficients block (after the count parameters)?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines, in this
#'   case to the \code{summary} method for the object.
#'
#' @method extract zeroinfl
#' @aliases extract.zeroinfl
#' @importFrom stats AIC logLik
#' @export
setMethod("extract", signature = className("zeroinfl", "pscl"),
          definition = extract.zeroinfl)


# -- extract.hurdle (pscl) -----------------------------------------------------

#' @noRd
extract.hurdle <- extract.zeroinfl

#' \code{\link{extract}} method for \code{hurdle} objects
#'
#' \code{\link{extract}} method for \code{hurdle} objects created by the
#' \code{hurdle} function in the \pkg{pscl} package.
#'
#' @inheritParams extract,zeroinfl-method
#'
#' @method extract hurdle
#' @aliases extract.hurdle
#' @importFrom stats AIC logLik
#' @export
setMethod("extract", signature = className("hurdle", "pscl"),
          definition = extract.hurdle)

# -- extract.logitr (logitr >= 0.8.0) ------------------------------------------

#' @noRd
extract.logitr <- function(model,
                           include.nobs = TRUE,
                           include.loglik = TRUE,
                           include.aic = TRUE,
                           include.bic = TRUE,
                           ...) {

  s <- summary(model, ...)

  coefnames <- names(s$coefficients)
  coefs <- s$coefTable$Estimate
  se <- s$coefTable$`Std. Error`
  pval <- s$coefTable$`Pr(>|z|)`

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.nobs == TRUE) {
    n <- s$n$obs
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.loglik == TRUE) {
    ll <- s$logLik
    gof <- c(gof, ll)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.aic == TRUE) {
    gof <- c(gof, s$statTable["AIC", ])
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    gof <- c(gof, s$statTable["BIC", ])
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

#' \code{\link{extract}} method for \code{logitr} objects
#'
#' \code{\link{extract}} method for \code{logitr} objects created by the
#' \code{logitr} function in the \pkg{logitr} package.
#'
#' @param model A statistical model object.
#' @param include.nobs Include the number of observations in summary table?
#' @param include.loglik Include the log-likelihood in summary table?
#' @param include.aic Include the the AIC in summary table?
#' @param include.bic Include the the BIC in summary table?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @method extract logitr
#' @aliases logitr
#' @author John Paul Helveston, \email{john.helveston@gmail.com}
#' @export
setMethod("extract", signature = className("logitr", "logitr"),
          definition = extract.logitr)
