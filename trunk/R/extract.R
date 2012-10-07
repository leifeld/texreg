# The texreg package was written by Philip Leifeld.
# Please use the forum at http://r-forge.r-project.org/projects/texreg/ 
# for bug reports, help or feature requests.


# generic extract function
setGeneric("extract", function(model, ...) standardGeneric("extract"), 
    package="texreg")


# extension for clogit objects (survival package); submitted by Sebastian Daza
extract.clogit <- function(model) {
  tab <- summary(model)$coef[,c(-2,-4)]
  aic <- extractAIC(model)[2]
  event <- model$nevent
  n <- model$n
  mis <- length(model$na.action)
  gof <- matrix(c(aic, event, n, mis), ncol=1)
  row.names(gof) <- c("AIC", "Events", "Num. obs.", "Missings")
  
  table.content <- list(tab, gof)
  return(table.content)
}

setMethod("extract", signature=className("clogit", "survival"), 
    definition = extract.clogit)


# extension for ergm objects
extract.ergm <- function(model) {
  tab <- summary(model)$coefs[,-3] #extract coefficient table
  
  lik <- model$mle.lik[1] #extract log likelihood
  aic <- summary(model)$aic #extract AIC
  bic <- summary(model)$bic #extract BIC
  gof <- matrix(c(aic, bic, lik), ncol=1)
  row.names(gof) <- c("AIC", "BIC", "Log Likelihood")
  
  table.content <- list(tab, gof)
  return(table.content)
}

setMethod("extract", signature=className("ergm", "ergm"), 
    definition = extract.ergm)


# extension for glm objects
extract.glm <- function(model) {
  tab <- summary(model)$coef[,-3] #extract coefficient table
  
  aic <- summary(model)$aic #extract AIC
  n <- nobs(model) #extract number of observations
  
  gof <- matrix(c(aic, n), ncol=1)
  row.names(gof) <- c("AIC", "Num. obs.")
  
  table.content <- list(tab, gof)
  return(table.content)
}

setMethod("extract", signature=className("glm", "stats"), 
    definition = extract.glm)


# extension for gls objects
extract.gls <- function(model) {
  tab <- summary(model)$tTable[,-3] #extract coefficient table
  
  lik <- summary(model)$logLik #extract log likelihood
  aic <- summary(model)$AIC #extract AIC
  bic <- summary(model)$BIC #extract BIC
  n <- nobs(model) #extract number of observations
  gof <- matrix(c(aic, bic, lik, n), ncol=1)
  row.names(gof) <- c("AIC", "BIC", "Log Likelihood", "Num. obs.")
  
  table.content <- list(tab, gof)
  return(table.content)
}

setMethod("extract", signature=className("gls", "nlme"), 
    definition = extract.gls)


# extension for lm objects
extract.lm <- function(model) {
  tab <- summary(model)$coef[,-3] #extract coefficient table
  
  rs <- summary(model)$r.squared #extract R-squared
  adj <- summary(model)$adj.r.squared #extract adjusted R-squared
  n <- nobs(model) #extract number of observations
  
  gof <- matrix(c(rs, adj, n), ncol=1)
  row.names(gof) <- c("R$^2$", "Adj. R$^2$", "Num. obs.")
  
  table.content <- list(tab, gof)
  return(table.content)
}

setMethod("extract", signature=className("lm", "stats"), 
    definition = extract.lm)


# extension for lme objects
extract.lme <- function(model) {
  tab <- summary(model)$tTable[,-3:-4] #extract coefficient table
  
  lik <- summary(model)$logLik #extract log likelihood
  aic <- summary(model)$AIC #extract AIC
  bic <- summary(model)$BIC #extract BIC
  n <- nobs(model) #extract number of observations
  gof <- matrix(c(aic, bic, lik, n), ncol=1)
  row.names(gof) <- c("AIC", "BIC", "Log Likelihood", "Num. obs.")
  
  table.content <- list(tab, gof)
  return(table.content)
}

setMethod("extract", signature=className("lme", "nlme"), 
    definition = extract.lme)


# extension for lnam objects (sna package)
extract.lnam <- function(model) {
  coefs <- coef(model)
  se <- se.lnam(fit)
  p <- 2 * (1 - pnorm(abs(coefs), 0, se))
  tab <- cbind(coefs, se, p)
  
  rss <- sum(model$residuals^2)
  mss <- sum((model$fitted - mean(model$fitted))^2)
  rdfns <- model$df.residual + 1
  rsquared <- mss / (mss + rss)
  adj.rsquared <- 1 - (1 - mss / (mss + rss)) * model$df.total / rdfns
  loglik <- model$lnlik.model
  aic <- -2 * model$lnlik.model + 2 * model$df.model
  bic <- -2 * model$lnlik.model + log(model$df.total) * model$df.model
  
  gof <- matrix(c(rsquared, adj.rsquared, loglik, aic, bic), ncol=1)
  row.names(gof) <- c("R$^2$", "Adj. R$^2$", "Log Likelihood", "AIC", "BIC")
  
  table.content <- list(tab, gof)
  return(table.content)
}

setMethod("extract", signature=className("lnam", "sna"), 
    definition = extract.lnam)


# extension for lrm objects (Design or rms package); submitted by Fabrice Le Lec
extract.lrm <- function(model) {
  tab <- model$coef #extract coefficient table
  
  attributes(model$coef)$names <- lapply(attributes(model$coef)$names, 
    function(x) gsub(">=", " $\\\\geq$ ", x))
  
  tab <- cbind(COEFEST = model$coef, SE = sqrt(diag(model$var)), 
      PVALUES = pnorm(abs(model$coef/sqrt(diag(model$var))), 
      lower.tail = FALSE)*2)
  
  pseudors <- model$stats[10] #extract pseudo R-squared
  LR <- model$stats[3] #extract LR
  n <- model$stats[1] #extract number of observations
  gof <- matrix(c(pseudors, LR, n), ncol=1)
  row.names(gof) <- c("Pseudo R$^2$", "L.R.", "Num. obs.")
  
  table.content <- list(tab, gof) #put coefficients and gofs in a list
  return(table.content) #return the list object
}

setMethod("extract", signature=className("lrm", "rms"), 
    definition = extract.lrm)
setMethod("extract", signature=className("lrm", "Design"), 
    definition = extract.lrm)


# extension for plm objects (from the plm package); submitted by Lena Koerber
extract.plm <- function(model) {
  tab <- summary(model)$coef[,-3]
  
  rs <- summary(model)$r.squared
  n <- length(summary(model)$resid)
  gof <- matrix(c(rs, n), ncol=1)
  row.names(gof) <- c("R$^2$", "Adj. R$^2$", "Num. obs.")
  
  table.content <- list(tab, gof)
  return(table.content)
}

setMethod("extract", signature=className("plm", "plm"), 
    definition = extract.plm)


# extension for pmg objects (from the plm package); submitted by Lena Koerber
extract.pmg <- function(model) {
  co <- data.matrix(summary(model)$coef)
  se <- (diag(summary(model)$vcov))^(1/2) #standard errors
  t <- co / se #t-statistics
  n <- length(summary(model)$resid) #number of observations
  d <- n - length(co) #degrees of freedom
  pval <- 2 * pt(-abs(t), df=d)
  tab <- cbind(co, se, pval) #coefficient table
  
  gof <- matrix(n, ncol=1)
  row.names(gof) <- c("Num. obs.")
  
  table.content <- list(tab, gof)
  return(table.content)
}

setMethod("extract", signature=className("pmg", "plm"), 
    definition = extract.pmg)


# extension for rq objects (quantreg package); submitted by Lena Koerber
extract.rq <- function(model) {
  tab <- summary(model, cov=TRUE)$coef[,-3]
  n <- length(summary(model)$resid)
  tau<-summary(model)$tau
  gof <- matrix(c(n, tau), ncol=1)
  row.names(gof) <- c("Num. obs.", "Percentile")
  
  table.content <- list(tab, gof)
  return(table.content)
}

setMethod("extract", signature=className("rq", "quantreg"), 
    definition = extract.rq)


# extension for systemfit objects; submitted by Johannes Kutsam
extract.systemfit <- function(model) {
  equationList <- list()
  for(eq in model$eq){  #go through estimated equations
    sum <- summary(eq)  #extract model summary
    tab <- coef(sum)[,-3]  #coefficients table for the equation
    
    rsquared <- sum$r.squared  #extract r-squared
    radj <- sum$adj.r.squared  #extract adjusted r-squared
    n <- nobs(model)  #extract number of observations
    
    gof <- matrix(c(rsquared, radj, n), ncol=1)  #GOF matrix
    row.names(gof) <- c("R$^2$", "Adj. R$^2$", "Num. obs.")  #set GOF row names

    table.content <- list(tab, gof)
    equationList[[eq$eqnNo]] <- table.content
  }
  return(equationList)  #returns a list of table.content lists
}

setMethod("extract", signature=className("systemfit", "systemfit"), 
    definition = extract.systemfit)

