## extension for feglm objects from alpaca package
## contributed by Oliver Reiter
## Github: zauster
## E: oliver_reiter@gmx.at
extract.feglm <- function(model,
                          include.nobs = TRUE,
                          include.deviance = TRUE,
                          ## include.loglik = TRUE,
                          ## include.aic = TRUE,
                          ## include.bic = TRUE,
                          ...) {
    s <- summary(model)
    names <- names(model$coefficients)
    co <- s$cm[, 1]
    se <- s$cm[, 2]
    pval <- s$cm[, 4]
    n <- s$nobs

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
        gof <- c(gof, n)
        gof.names <- c(gof.names, "Num. obs.")
        gof.decimal <- c(gof.decimal, FALSE)
    }
    ## if (include.loglik == TRUE) {
    ## }
    ## if (include.aic == TRUE) {
    ## }
    ## if (include.bic == TRUE) {
    ## }
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

setMethod("extract", signature = className("feglm", "alpaca"),
          definition = extract.feglm)
