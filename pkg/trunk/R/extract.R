# The texreg package was written by Philip Leifeld.
# Please use the forum at http://r-forge.r-project.org/projects/texreg/ 
# for bug reports, help or feature requests.


# generic extract function
setGeneric("extract", function(model, ...) standardGeneric("extract"), 
    package="texreg")


# extension for clogit objects (survival package); submitted by Sebastian Daza
extract.clogit <- function(model) {
  coefficient.names <- rownames(summary(model)$coef)
  coefficients <- summary(model)$coef[,1]
  standard.errors <- summary(model)$coef[,3]
  significance <- summary(model)$coef[,5]
  
  aic <- extractAIC(model)[2]
  event <- model$nevent
  n <- model$n
  mis <- length(model$na.action)
  gof <- c(aic, event, n, mis)
  gof.names <- c("AIC", "Events", "Num. obs.", "Missings")
  
  tr <- createTexreg(
      coef.names=coefficient.names, 
      coef=coefficients, 
      se=standard.errors, 
      pvalues=significance, 
      gof.names=gof.names, 
      gof=gof
  )
  return(tr)
}

setMethod("extract", signature=className("clogit", "survival"), 
    definition = extract.clogit)


# extension for ergm objects
extract.ergm <- function(model) {
  coefficient.names <- rownames(summary(model)$coefs)
  coefficients <- summary(model)$coefs[,1]
  standard.errors <- summary(model)$coefs[,2]
  significance <- summary(model)$coefs[,4]
  
  lik <- model$mle.lik[1] #extract log likelihood
  aic <- summary(model)$aic #extract AIC
  bic <- summary(model)$bic #extract BIC
  gof <- c(aic, bic, lik)
  gof.names <- c("AIC", "BIC", "Log Likelihood")
  
  tr <- createTexreg(
      coef.names=coefficient.names, 
      coef=coefficients, 
      se=standard.errors, 
      pvalues=significance, 
      gof.names=gof.names, 
      gof=gof
  )
  return(tr)
}

setMethod("extract", signature=className("ergm", "ergm"), 
    definition = extract.ergm)


# extension for glm objects
extract.glm <- function(model) {
  coefficient.names <- rownames(summary(model)$coef)
  coefficients <- summary(model)$coef[,1]
  standard.errors <- summary(model)$coef[,2]
  significance <- summary(model)$coef[,4]
  
  aic <- summary(model)$aic #extract AIC
  n <- nobs(model) #extract number of observations
  
  gof <- c(aic, n)
  gof.names <- c("AIC", "Num. obs.")
  
  tr <- createTexreg(
      coef.names=coefficient.names, 
      coef=coefficients, 
      se=standard.errors, 
      pvalues=significance, 
      gof.names=gof.names, 
      gof=gof
  )
  return(tr)
}

setMethod("extract", signature=className("glm", "stats"), 
    definition = extract.glm)


# extension for gls objects
extract.gls <- function(model) {
  coefficient.names <- rownames(summary(model)$tTable)
  coefficients <- summary(model)$tTable[,1]
  standard.errors <- summary(model)$tTable[,2]
  significance <- summary(model)$tTable[,4]
  
  lik <- summary(model)$logLik #extract log likelihood
  aic <- summary(model)$AIC #extract AIC
  bic <- summary(model)$BIC #extract BIC
  n <- nobs(model) #extract number of observations
  gof <- c(aic, bic, lik, n)
  gof.names <- c("AIC", "BIC", "Log Likelihood", "Num. obs.")
  
  tr <- createTexreg(
      coef.names=coefficient.names, 
      coef=coefficients, 
      se=standard.errors, 
      pvalues=significance, 
      gof.names=gof.names, 
      gof=gof
  )
  return(tr)
}

setMethod("extract", signature=className("gls", "nlme"), 
    definition = extract.gls)


# extension for lm objects
extract.lm <- function(model) {
  names <- rownames(summary(model)$coef)
  co <- summary(model)$coef[,1]
  se <- summary(model)$coef[,2]
  pval <- summary(model)$coef[,4]
  
  rs <- summary(model)$r.squared #extract R-squared
  adj <- summary(model)$adj.r.squared #extract adjusted R-squared
  n <- nobs(model) #extract number of observations
  
  gof <- c(rs, adj, n)
  gof.names <- c("R$^2$", "Adj. R$^2$", "Num. obs.")
  
  tr <- createTexreg(
      coef.names=names, 
      coef=co, 
      se=se, 
      pvalues=pval, 
      gof.names=gof.names, 
      gof=gof
  )
  return(tr)
}

setMethod("extract", signature=className("lm", "stats"), 
    definition = extract.lm)


# extension for lme objects
extract.lme <- function(model) {
  coefficient.names <- rownames(summary(model)$tTable)
  coefficients <- summary(model)$tTable[,1]
  standard.errors <- summary(model)$tTable[,2]
  significance <- summary(model)$tTable[,5]
  
  lik <- summary(model)$logLik #extract log likelihood
  aic <- summary(model)$AIC #extract AIC
  bic <- summary(model)$BIC #extract BIC
  n <- nobs(model) #extract number of observations
  gof <- c(aic, bic, lik, n)
  gof.names <- c("AIC", "BIC", "Log Likelihood", "Num. obs.")
  
  tr <- createTexreg(
      coef.names=coefficient.names, 
      coef=coefficients, 
      se=standard.errors, 
      pvalues=significance, 
      gof.names=gof.names, 
      gof=gof
  )
  return(tr)
}

setMethod("extract", signature=className("lme", "nlme"), 
    definition = extract.lme)


# extension for lnam objects (sna package)
extract.lnam <- function(model) {
  coefs <- coef(model)
  coef.names <- names(coefs)
  se <- c(model$beta.se, model$rho1.se, model$rho2.se)
  p <- 2 * (1 - pnorm(abs(coefs), 0, se))
  
  rss <- sum(model$residuals^2)
  mss <- sum((model$fitted - mean(model$fitted))^2)
  rdfns <- model$df.residual + 1
  rsquared <- mss / (mss + rss)
  adj.rsquared <- 1 - (1 - mss / (mss + rss)) * model$df.total / rdfns
  loglik <- model$lnlik.model
  aic <- -2 * model$lnlik.model + 2 * model$df.model
  bic <- -2 * model$lnlik.model + log(model$df.total) * model$df.model
  
  gof <- c(rsquared, adj.rsquared, loglik, aic, bic)
  gof.names <- c("R$^2$", "Adj. R$^2$", "Log Likelihood", "AIC", "BIC")
  
  tr <- createTexreg(
      coef.names=coef.names, 
      coef=coefs, 
      se=se, 
      pvalues=p, 
      gof.names=gof.names, 
      gof=gof
  )
  return(tr)
}

setMethod("extract", signature=className("lnam", "sna"), 
    definition = extract.lnam)


# extension for lrm objects (Design or rms package); submitted by Fabrice Le Lec
extract.lrm <- function(model) {
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
  gof <- c(pseudors, LR, n)
  gof.names <- c("Pseudo R$^2$", "L.R.", "Num. obs.")
  
  tr <- createTexreg(
      coef.names=coef.names, 
      coef=coef, 
      se=se, 
      pvalues=p, 
      gof.names=gof.names, 
      gof=gof
  )
  return(tr)
}

setMethod("extract", signature=className("lrm", "rms"), 
    definition = extract.lrm)
setMethod("extract", signature=className("lrm", "Design"), 
    definition = extract.lrm)


# extension for plm objects (from the plm package); submitted by Lena Koerber
extract.plm <- function(model) {
  coefficient.names <- rownames(summary(model)$coef)
  coefficients <- summary(model)$coef[,1]
  standard.errors <- summary(model)$coef[,2]
  significance <- summary(model)$coef[,4]
  
  rs <- summary(model)$r.squared
  n <- length(summary(model)$resid)
  gof <- c(rs, n)
  gof.names <- c("R$^2$", "Adj. R$^2$", "Num. obs.")
  
  tr <- createTexreg(
      coef.names=coefficient.names, 
      coef=coefficients, 
      se=standard.errors, 
      pvalues=significance, 
      gof.names=gof.names, 
      gof=gof
  )
  return(tr)
}

setMethod("extract", signature=className("plm", "plm"), 
    definition = extract.plm)


# extension for pmg objects (from the plm package); submitted by Lena Koerber
extract.pmg <- function(model) {
  co <- summary(model)$coef
  se <- (diag(summary(model)$vcov))^(1/2) #standard errors
  t <- co / se #t-statistics
  n <- length(summary(model)$resid) #number of observations
  d <- n - length(co) #degrees of freedom
  pval <- 2 * pt(-abs(t), df=d)
  tab <- cbind(co, se, pval) #coefficient table
  names <- rownames(tab)
  
  gof <- n
  gof.names <- "Num. obs."
  
  tr <- createTexreg(
      coef.names=names, 
      coef=co, 
      se=se, 
      pvalues=pval, 
      gof.names=gof.names, 
      gof=gof
  )
  return(tr)
}

setMethod("extract", signature=className("pmg", "plm"), 
    definition = extract.pmg)


# extension for rq objects (quantreg package); submitted by Lena Koerber
extract.rq <- function(model) {
  co <- summary(model, cov=TRUE)$coef[,1]
  names <- rownames(summary(model, cov=TRUE)$coef)
  se <- summary(model, cov=TRUE)$coef[,2]
  pval <- summary(model, cov=TRUE)$coef[,4]
  
  n <- length(summary(model)$resid)
  tau<-summary(model)$tau
  gof <- c(n, tau)
  gof.names <- c("Num. obs.", "Percentile")
  
  tr <- createTexreg(
      coef.names=names, 
      coef=co, 
      se=se, 
      pvalues=pval, 
      gof.names=gof.names, 
      gof=gof
  )
  return(tr)
}

setMethod("extract", signature=className("rq", "quantreg"), 
    definition = extract.rq)


# extension for systemfit objects; submitted by Johannes Kutsam
extract.systemfit <- function(model) {
  equationList <- list()
  for(eq in model$eq){  #go through estimated equations
    sum <- summary(eq)  #extract model summary
    names <- rownames(coef(sum))
    co <- coef(sum)[,1]
    se <- coef(sum)[,2]
    pval <- coef(sum)[,4]
    
    rsquared <- sum$r.squared  #extract r-squared
    radj <- sum$adj.r.squared  #extract adjusted r-squared
    n <- nobs(model)  #extract number of observations
    
    gof <- c(rsquared, radj, n)
    gof.names <- c("R$^2$", "Adj. R$^2$", "Num. obs.")

    tr <- createTexreg(
      coef.names=names, 
      coef=co, 
      se=se, 
      pvalues=pval, 
      gof.names=gof.names, 
      gof=gof
    )
    equationList[[eq$eqnNo]] <- tr
  }
  return(equationList)  #returns a list of table.content lists
}

setMethod("extract", signature=className("systemfit", "systemfit"), 
    definition = extract.systemfit)

