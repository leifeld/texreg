# The texreg package was written by Philip Leifeld.
# Please use the forum at http://r-forge.r-project.org/projects/texreg/ 
# for bug reports, help or feature requests.


# a class for texreg objects
setClass(Class="texreg", 
    representation=representation(
        coef.names="character",  #row names of the coefficient block
        coef="numeric",          #the coefficients
        se="numeric",            #standard errors
        pvalues="numeric",       #p-values
        gof.names="character",   #row names of the goodness-of-fit block
        gof="numeric",           #goodness-of-fit statistics
        gof.decimal="logical"    #number of decimal places for each GOF value
    ), 
    validity=function(object) {
        if (length(object@coef.names) != length(object@coef) || 
            length(object@coef.names) != length(object@se) ||
            length(object@coef) != length(object@se)
        ) {
            stop("coef.names, coef, and se must have the same length!")
        }
        if (length(object@pvalues) != 0 && 
            (length(object@coef.names) != length(object@pvalues) || 
            length(object@coef) != length(object@pvalues) || 
            length(object@se) != length(object@pvalues))
        ) {
          stop(paste("pvalues must have the same length as coef.names, coef,", 
              "and se, or it must have a length of zero."))
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
        return(TRUE)
    }
)

# constructor for texreg objects
createTexreg <- function(coef.names, coef, se, pvalues=numeric(0), 
    gof.names=character(0), gof=numeric(0), gof.decimal=logical(0)) {
  new("texreg", coef.names=coef.names, coef=coef, se=se, pvalues=pvalues, 
      gof.names=gof.names, gof=gof, gof.decimal=gof.decimal)
}

# define show method for pretty output of texreg objects
setMethod(f="show", signature="texreg", definition=function(object) {
  if (length(object@pvalues) > 0) {
    cat("\n")
    coefBlock <- cbind(object@coef, object@se, object@pvalues)
    colnames(coefBlock) <- c("coef.", "s.e.", "p")
  } else {
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
    gofBlock <- matrix(object@gof, ncol=1)
    colnames(gofBlock) <- "GOF"
  } else {
    gofBlock <- data.frame(cbind(object@gof,object@gof.decimal), ncol=2)
    colnames(gofBlock) <- c("GOF", "dec. places")
  }
  rownames(gofBlock) <- object@gof.names
  print(gofBlock)
  cat("\n")
})

