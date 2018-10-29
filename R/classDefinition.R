# The texreg package was written by Philip Leifeld.
# Please use the forum at http://r-forge.r-project.org/projects/texreg/ 
# for bug reports, help or feature requests.


# a class for texreg objects
setClass(Class = "texreg", 
    representation = representation(
        coef.names = "character",  # row names of the coefficient block
        coef = "numeric",          # the coefficients
        se = "numeric",            # standard errors
        pvalues = "numeric",       # p-values
        ci.low = "numeric",        # lower confidence interval
        ci.up = "numeric",         # upper confidence interval
        gof.names = "character",   # row names of the goodness-of-fit block
        gof = "numeric",           # goodness-of-fit statistics
        gof.decimal = "logical",   # number of decimal places for each GOF value
        model.name = "character",  # name of the model
        Df = "numeric"           # residuals degrees of freedom
    ), 
    validity = function(object) {
        if (length(object@coef.names) != length(object@coef)) {
            stop("coef.names and coef must have the same length!")
        }
        if (length(object@pvalues) != 0 && 
            (length(object@coef.names) != length(object@pvalues) || 
            length(object@coef) != length(object@pvalues))
        ) {
          stop(paste("pvalues must have the same length as coef.names, coef,", 
              "and se, or it must have a length of zero."))
        }
        if ((length(object@ci.low) != length(object@ci.up)) || 
            (length(object@ci.low) != 0 && length(object@ci.low) != 
            length(object@coef))) {
          stop("CIs must have a length of zero or the same length as coef.")
        }
        if (length(object@gof.names) != length(object@gof)) {
            stop("gof.names and gof must have the same length!")
        }
        if (length(object@gof.decimal) != 0 && 
            (length(object@gof.decimal) != length(object@gof) || 
            length(object@gof.decimal) != length(object@gof.names))
        ) {
          stop(paste("gof.decimal must have the same length as gof.names and", 
              "gof, or it must have a length of zero."))
        }
        if (length(object@model.name) > 1) {
          stop("Only one model name can be provided.")
        }
        return(TRUE)
    }
)

# constructor for texreg objects
createTexreg <- function(coef.names, coef, se = numeric(0), 
    pvalues = numeric(0), ci.low = numeric(0), ci.up = numeric(0), 
    gof.names = character(0), gof = numeric(0), gof.decimal = logical(0), 
    model.name = character(0), Df = numeric(0)) {
  new("texreg", coef.names = coef.names, coef = coef, se = se, 
      pvalues = pvalues, ci.low = ci.low, ci.up = ci.up, gof.names = gof.names, 
      gof = gof, gof.decimal = gof.decimal, model.name = model.name, Df = Df)
}

# define show method for pretty output of texreg objects
setMethod(f = "show", signature = "texreg", definition = function(object) {
  if (length(object@model.name) == 1) {
    cat(paste("Model name:", object@model.name))
  }
  if (length(object@se) == 0 && length(object@ci.up) > 0) {
    coefBlock <- cbind(object@coef, object@ci.low, object@ci.up)
    colnames(coefBlock) <- c("coef.", "lower CI", "upper CI")
  } else if (length(object@se) == 0 && length(object@pvalues) == 0) {
    cat(paste("\nNo standard errors and p-values were defined for this", 
          "texreg object.\n"))
    coefBlock <- cbind(object@coef)
    colnames(coefBlock) <- "coef."
  } else if (length(object@se) == 0) {
    cat(paste("\nNo standard errors were defined for this texreg object.\n"))
    coefBlock <- cbind(object@coef, object@pvalues)
    colnames(coefBlock) <- c("coef.", "p")
  } else if (length(object@pvalues) > 0) {
    cat("\n")
    coefBlock <- cbind(object@coef, object@se, object@pvalues)
    colnames(coefBlock) <- c("coef.", "s.e.", "p")
  } else if (length(object@Df)  == 0) {
      cat(paste("\nNo degrees of freedom were defined for this texreg object.\n"))
  } else if (length(object@Df)  == 1) {
      coefBlock <- cbind(object@Df)
      colnames(coefBlock) <- c("DF")
  }else {
    cat("\nNo p-values were defined for this texreg object.\n")
    coefBlock <- cbind(object@coef, object@se)
    colnames(coefBlock) <- c("coef.", "s.e.")
  }
  rownames(coefBlock) <- object@coef.names
  dec <- object@gof.decimal
  if (length(dec) == 0) {
    cat("No decimal places were defined for the GOF statistics.\n")
  }
  cat("\n")
  print(coefBlock)
  cat("\n")
  if (length(dec) == 0) {
    gofBlock <- matrix(object@gof, ncol = 1)
    colnames(gofBlock) <- "GOF"
  } else {
    gofBlock <- data.frame(object@gof,object@gof.decimal)
    colnames(gofBlock) <- c("GOF", "dec. places")
  }
  rownames(gofBlock) <- object@gof.names
  if (nrow(gofBlock) > 0) {
    print(gofBlock)
    cat("\n")
  } else {
    cat("No GOF block defined.\n")
  }
})