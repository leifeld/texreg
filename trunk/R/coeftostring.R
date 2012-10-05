# The texreg package was written by Philip Leifeld.
# Please use the forum at http://r-forge.r-project.org/projects/texreg/ 
# for bug reports, help or feature requests.


# function which reformats a coefficient with two decimal places
coeftostring <- function(x, lead.zero=FALSE, digits=2) {
  if (!is.finite(x)) {
    return("")
  }
  if (digits < 0) {
    stop("The number of digits must be 0 or higher.")
  }
  y <- as.character(round(x, digits))
  
  actual.digits <- attributes(regexpr("\\.[0-9]+$", y))$match.length - 1
  if (grepl("\\.", y) == FALSE) { #no decimal places present
    zeros <- paste(rep("0", digits), collapse="")
    y <- paste(y, ".", zeros, sep="")
  } else if (grepl("\\.[0-9]$", y) == TRUE && digits > 1) { #only one decimal p.
    zeros <- paste(rep("0", digits-1), collapse="")
    y <- paste(y, zeros, sep="")
  } else if (actual.digits < digits) { #more desired digits than present
    fill <- digits - actual.digits
    zeros <- paste(rep("0", fill), collapse="")
    y <- paste(y, zeros, sep="")
  }
  if (grepl("\\.$", y) == TRUE) {
    y <- gsub("\\.$", "", y)
  }
  if (lead.zero==FALSE && (grepl("^0",y) == TRUE || grepl("^-0",y) == TRUE)) {
    y <- gsub("0\\.", "\\.", y)
  }
  return(y)
}

