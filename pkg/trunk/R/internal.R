# The texreg package was written by Philip Leifeld.
# Please use the forum at http://r-forge.r-project.org/projects/texreg/ 
# for bug reports, help or feature requests.


# display version number and date when the package is loaded
.onAttach <- function(libname, pkgname) {
  desc  <- packageDescription(pkgname, libname)
  packageStartupMessage('Version:  ', desc$Version, '\n', 'Date:     ', 
      desc$Date)
}


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


# function which conflates a matrix with duplicate row names
rearrangeMatrix <- function(m) {
  
  # The following code block rearranges a matrix with duplicate row names such 
  # that these rows are conflated where possible. First, an empty matrix q with
  # the same width is created. The rows will be copied iteratively into this 
  # matrix. Second, we go through the unique row names, and for each row name 
  # we create a small virtual matrix in which the values will be nicely 
  # rearranged. After rearranging the values, this small matrix is rbinded to 
  # the q matrix. Rearranging works in the following way (the inner loop): for 
  # every column, we create a vector of all values corresponding to the specific
  # row name (as specified by the outer loop). We retain only non-NA values 
  # because irrelevant information should be removed from the coefficients 
  # table. Then we put the first non-NA value in the first vertical slot of the 
  # virtual matrix, the second non-NA value of the same row name in the second 
  # slot, etc., and we create additional rows in the virtual matrix as needed.
  # By doing this, we ensure that no space in the matrix is wasted with NA 
  # values. When going to the next column, we place the non-NA values in the 
  # correct slot again, and we only create new rows if needed. The virtual rows 
  # are finally rbinded to the large replacement matrix q.
  
  unique.names <- unique(rownames(m))              #unique row names in m
  num.unique <- length(unique.names)               #count these unique names
  orig.width <- length(m[1,])                      #number of columns in m
  q <- matrix(nrow=0, ncol=orig.width)             #new matrix with same width
  for (i in 1:num.unique) {                        #go through unique row names
    rows <- matrix(NA, nrow=0, ncol=orig.width)    #create matrix where re-
                                                   #arranged rows will be stored
    for (j in 1:orig.width) {                      #go through columns in m
      current.name <- unique.names[i]              #save row name
      nonNa <- m[rownames(m)==current.name,j]      #create a vector of values
                                                   #with same rowname in the col
      nonNa <- nonNa[!is.na(nonNa)]                #retain only non-NA values
      for (k in 1:length(nonNa)) {                 #go through non-NA values
        if (k > dim(rows)[1]) {                    #add an NA-only row in which
          rows <- rbind(rows, rep(NA, orig.width)) #the values are stored
          rownames(rows)[k] <- unique.names[i]     #also add the row name
        }
        rows[k,j] <- nonNa[k]                      #actually store the value
      }
    }
    q <- rbind(q, rows)                            #add the new row(s) to q
  }
  return(q)
}


# function which wraps models in a list and extracts texreg objects from them
get.data <- function(l, ...) {

  # if a single model is handed over, put model inside a list
  if (!"list" %in% class(l)) {
    l <- list(l)
  }

  # extract data from the models
  models <- NULL
  for (i in 1:length(l)) {
    model <- extract(l[[i]], ...)
    if (class(model) == "list") {       #nested list of models (e.g. systemfit)
      models <- append(models, model)
    } else {                            #normal case; one model
      models <- append(models, list(model))
    }
  }
  
  return(models)
}


# function which extracts names of the goodness-of-fit statistics
get.gof <- function(models) {
  gof.names <- character()  #names of all models in one vector
  for (i in 1:length(models)) {
    gn <- models[[i]]@gof.names
    for (j in 1:length(gn)) {
      if (!gn[j] %in% gof.names) {
        gof.names <- append(gof.names, gn[j])
      }
    }
  }
  return(gof.names)
}


# function which replaces coefs, SEs and p values by custom values if provided
override <- function(models, override.coef, override.se, override.pval) {
  
  for (i in 1:length(models)) {
    
    # coefficients
    if (class(override.coef) != "list" && length(override.coef) == 1 && 
        override.coef == 0) {
      cf <- models[[i]]@coef
    } else if (class(override.coef) == "numeric" && length(models) == 1 && 
        length(override.coef) == length(models[[i]]@coef)) {
      cf <- override.coef
    } else if (class(override.coef) != "list") {
      warning("Coefficients must be provided as a list. Using default values.")
      cf <- models[[i]]@coef
    } else if (length(override.coef) != length(models)) {
      warning(paste("Number of coefficients provided does not match number of", 
          "models. Using default values."))
      cf <- models[[i]]@coef
    } else if (length(models[[i]]@coef) != length(override.coef[[i]])) {
      warning(paste("Number of coefficients provided does not match number of ",
          "terms in model ", i, ". Using default values.", sep=""))
      cf <- models[[i]]@coef
    } else if (class(override.coef[[i]]) != "numeric") {
      warning(paste("Coefficients provided for model", i, 
          "are not numeric. Using default values."))
      cf <- models[[i]]@coef
    } else {
      cf <- override.coef[[i]]
    }
    models[[i]]@coef <- cf
    
    # standard errors
    if (class(override.se) != "list" && length(override.se) == 1 && 
        override.se == 0) {
      se <- models[[i]]@se
    } else if (class(override.se) == "numeric" && length(models) == 1 && 
        length(override.se) == length(models[[i]]@se)) {
      se <- override.se
    } else if (class(override.se) != "list") {
      warning("SEs must be provided as a list. Using default SEs.")
      se <- models[[i]]@se
    } else if (length(override.se) != length(models)) {
      warning(paste("Number of SEs provided does not match number of models.", 
          "Using default SEs."))
      se <- models[[i]]@se
    } else if (length(models[[i]]@se) != length(override.se[[i]])) {
      warning(paste("Number of SEs provided does not match number of ", 
          "coefficients in model ", i, ". Using default SEs.", sep=""))
      se <- models[[i]]@se
    } else if (class(override.se[[i]]) != "numeric") {
      warning(paste("SEs provided for model", i, 
          "are not numeric. Using default SEs."))
      se <- models[[i]]@se
    } else {
      se <- override.se[[i]]
    }
    models[[i]]@se <- se
    
    # p values
    if (class(override.pval) != "list" && length(override.pval) == 1 && 
        override.pval == 0) {
      pval <- models[[i]]@pvalues
    } else if (class(override.pval) == "numeric" && length(models) == 1 && 
        length(override.pval) == length(models[[i]]@pvalues)) {
      pval <- override.pval
    } else if (class(override.pval) != "list") {
      warning("p values must be provided as a list. Using default p values.")
      pval <- models[[i]]@pvalues
    } else if (length(override.pval) != length(models)) {
      warning(paste("Number of p values provided does not match number of", 
          "models. Using default p values."))
      pval <- models[[i]]@pvalues
    } else if (length(models[[i]]@se) != length(override.pval[[i]])) {
      # previous line: comparison with se because pvalues can be empty
      warning(paste("Number of p values provided does not match number of ", 
          "coefficients in model ", i, ". Using default p values.", sep=""))
      pval <- models[[i]]@pvalues
    } else if (class(override.pval[[i]]) != "numeric") {
      warning(paste("p values provided for model", i, 
          "are not numeric. Using default p values."))
      pval <- models[[i]]@pvalues
    } else {
      pval <- override.pval[[i]]
    }
    models[[i]]@pvalues <- pval
  }
  
  return(models)
}


# function which converts LaTeX code in GOF names to HTML oder text/screen code
tex.replace <- function(models, type="html") {
  for (i in 1:length(models)) {
    if (type == "html") {
      r <- "<sup>2</sup>"
    } else if (type == "screen") {
      r <- "^2"
    }
    models[[i]]@gof.names <- sub("\\$\\^2\\$", r, models[[i]]@gof.names)
    models[[i]]@gof.names <- sub("\\\\ ", " ", models[[i]]@gof.names)
    models[[i]]@gof.names <- sub("\\ ", " ", models[[i]]@gof.names)
  }
  return(models)
}


# put models and GOFs into a common matrix
aggregate.matrix <- function(models, gof.names, custom.gof.names, digits, 
    returnobject="m") {

  # aggregate GOF statistics in a matrix and create list of coef blocks
  gofs <- matrix(nrow=length(gof.names), ncol=length(models))
  row.names(gofs) <- gof.names
  coefs <- list()
  decimal.matrix <- matrix(nrow=length(gof.names), ncol=length(models))
  for (i in 1:length(models)) {
    cf <- models[[i]]@coef
    se <- models[[i]]@se
    pv <- models[[i]]@pvalues
    if (length(pv) > 0) {
      coef <- cbind(cf, se, pv)
    } else { #p-values not provided -> use p-values of 0.99
      coef <- cbind(cf, se, rep(0.99, length(cf)))
    }
    rownames(coef) <- models[[i]]@coef.names
    coefs[[i]] <- coef
    for (j in 1:length(models[[i]]@gof)) {
      rn <- models[[i]]@gof.names[j]
      val <- models[[i]]@gof[j]
      col <- i
      if (is.na(models[[i]]@gof.decimal[j])) {
        dec <- digits
      } else if (models[[i]]@gof.decimal[j] == FALSE) {
        dec <- 0
      } else {
        dec <- digits
      }
      row <- which(row.names(gofs) == rn)
      gofs[row,col] <- val
      decimal.matrix[row,col] <- dec
    }
  }
  
  # figure out correct order of the coefficients
  coef.order <- character()
  for (i in 1:length(coefs)) {
    for (j in 1:length(rownames(coefs[[i]]))) {
      if (!rownames(coefs[[i]])[j] %in% coef.order) {
        coef.order <- append(coef.order, rownames(coefs[[i]])[j])
      }
    }
  }
  
  # merge the coefficient tables
  if (length(coefs) == 1) {
    m <- coefs[[1]]
  } else if (length(coefs) > 1) {
    m <- coefs[[1]]
    for (i in 2:length(coefs)) {
      m <- merge(m, coefs[[i]], by=0, all=TRUE)
      rownames(m) <- m[,1]
      m <- m[,colnames(m)!="Row.names"]
      colnames(m) <- NULL
    }
  }
  colnames(m) <- rep(colnames(coefs[[1]]), length(coefs))
  
  # reorder merged coefficient table
  m.temp <- matrix(nrow=nrow(m), ncol=ncol(m))
  for (i in 1:nrow(m)) {
    new.row <- which(coef.order == rownames(m)[i])
    for (j in 1:length(m[i,])) {
      m.temp[new.row,j] <- m[i,j]
    }
  }
  rownames(m.temp) <- coef.order
  colnames(m.temp) <- colnames(m)
  m <- m.temp
  
  if (returnobject == "m") {
    return(m)
  } else if (returnobject == "gofs") {
  
    #replace GOF names by custom names
    if (is.null(custom.gof.names)) {
      #do nothing
    } else if (class(custom.gof.names) != "character") {
      stop("Custom GOF names must be provided as a vector of strings.")
    } else if (length(custom.gof.names) != length(gof.names)) {
      stop(paste("There are", length(gof.names), 
          "GOF statistics, but you provided", length(custom.gof.names), 
          "custom names for them."))
    } else if (any(is.na(custom.gof.names))) {
      stop("Custom GOF names are not allowed to contain NA values.")
    } else {
      rownames(gofs) <- custom.gof.names
    }
    
    return(gofs)
    
  } else if (returnobject == "decimal.matrix") {
    return(decimal.matrix)
  }
}


# use custom coefficient names if provided
customnames <- function(m, custom.names) {
  if (is.null(custom.names)) {
    return(m)
  } else if (length(custom.names) > 1) {
    if (!class(custom.names) == "character") {
      stop("Custom coefficient names must be provided as a vector of strings!")
    } else if (length(custom.names) != length(rownames(m))) {
      stop(paste("There are", length(rownames(m)), 
          "coefficients, but you provided", length(custom.names), 
          "custom names for them."))
    } else {
      rownames(m) <- custom.names
    }
  } else if (!is.na(custom.names) & class(custom.names) != "character") {
    stop("Custom coefficient names must be provided as a vector of strings.")
  } else if (length(custom.names) == 1 & class(custom.names) == "character") {
    rownames(m) <- custom.names
  }
  return(m)
}

## use custom GOF names if provided
#customgof <- function(gof.names, custom.gof.names) {
#  if (is.null(custom.gof.names)) {
#    return(gof.names)
#  } else if (class(custom.gof.names) != "character") {
#    stop("Custom GOF names must be provided as a vector of strings.")
#  } else if (length(custom.gof.names) != length(gof.names)) {
#    stop(paste("There are", length(gof.names), 
#        "GOF statistics, but you provided", length(custom.gof.names), 
#        "custom names for them."))
#  } else if (any(is.na(custom.gof.names))) {
#    stop("Custom GOF names are not allowed to contain NA values.")
#  } else {
#    return(custom.gof.names)
#  }
#  return(m)
#}


# remove coefficient rows that match the omit.coef regular expression
omitcoef <- function(m, omit.coef) {
  if (!is.na(omit.coef)) {
    if (!is.character(omit.coef)) {
      stop("omit.coef must be a character string!")
    }
    remove.rows <- grep(omit.coef, rownames(m))
    m <- m[-remove.rows,]
  }
  return(m)
}


# decide if default or custom model names should be used and return them
modelnames <- function(models, model.names) {
  if (is.null(model.names)) {
    return(paste("Model", 1:length(models)))
  } else if (length(model.names) > 1) {
    if (class(model.names) != "character") {
      stop("Model names must be specified as a vector of strings.")
    } else if (length(model.names) != length(models)) {
      stop(paste("There are", length(models), "models, but you provided", 
          length(model.names), "names for them."))
    } else {
      return(model.names)
    }
  } else if (!is.na(model.names) & class(model.names) != "character") {
    stop("Model names must be specified as a vector of strings.")
  } else if (class(model.names) == "character" & 
      length(model.names) != length(models)) {
    stop(paste("A single model name was specified. But there are in fact", 
        length(models), "models."))
  } else if (class(model.names) == "character") {
    return(model.names)
  } else {
    return(paste("Model", 1:length(models)))
  }
}


# return the output matrix with coefficients, SEs and significance stars
outputmatrix <- function(m, single.row, neginfstring, leading.zero, digits, 
    se.prefix, se.suffix, star.prefix, star.suffix, strong.signif, 
    stars, dcolumn=TRUE, symbol) {

  # write coefficient rows
  if (single.row==TRUE) {
    output.matrix <- matrix(ncol=(length(m)/3)+1, nrow=length(m[,1]))
    
    # row labels
    for (i in 1:length(rownames(m))) {
      output.matrix[i,1] <- rownames(m)[i]
    }
    
    # coefficients and standard errors
    for (i in 1:length(m[,1])) { #go through rows
      j <- 1 #column in the original, merged coef table
      k <- 2 #second column of output.matrix, i.e., coefficients
      while (j <= length(m)) {
        if (is.na(m[i,j])) {
          output.matrix[i,k] <- ""
        } else if (m[i,j] == -Inf) {
          output.matrix[i,k] <- neginfstring
        } else {
          std <- paste(se.prefix, coeftostring(m[i,j+1], leading.zero, 
              digits=digits), se.suffix, sep="")
          if (strong.signif == TRUE && stars==TRUE) {
            if (m[i,j+2] <= 0.001) {
              p <- paste(star.prefix, "***", star.suffix, sep="")
            } else if (m[i,j+2] <= 0.01) {
              p <- paste(star.prefix, "**", star.suffix, sep="")
            } else if (m[i,j+2] <= 0.05) {
              p <- paste(star.prefix, "*", star.suffix, sep="")
            } else if (m[i,j+2] <= 0.1) {
              p <- paste(star.prefix, symbol, star.suffix, sep="")
            } else {
              p <- ""
            }
          } else if (stars==TRUE) {
            if (m[i,j+2] <= 0.01) {
              p <- paste(star.prefix, "***", star.suffix, sep="")
            } else if (m[i,j+2] <= 0.05) {
              p <- paste(star.prefix, "**", star.suffix, sep="")
            } else if (m[i,j+2] <= 0.1) {
              p <- paste(star.prefix, "*", star.suffix, sep="")
            } else {
              p <- ""
            }
          } else {
              p <- ""
          }
          if (dcolumn == TRUE) {
            dollar <- ""
          } else {
            dollar <- "$"
          }
          entry <- paste(dollar, coeftostring(m[i,j], leading.zero, 
              digits=digits), std, p, dollar, sep="")
          output.matrix[i,k] <- entry
          
        }
        k <- k+1
        j <- j+3
      }
    }
  } else {
    output.matrix <- matrix(ncol=(length(m)/3)+1, nrow=2*length(m[,1]))
    
    # row labels
    for (i in 1:length(rownames(m))) {
      output.matrix[(i*2)-1,1] <- rownames(m)[i]
      output.matrix[(i*2),1] <- ""
    }
    
    # coefficients and standard deviations
    for (i in 1:length(m[,1])) {
      j <- 1
      k <- 2
      while (j <= length(m)) {
        if (is.na(m[i,j]) || is.nan(m[i,j])) {
          output.matrix[(i*2)-1,k] <- "" #upper coefficient row
          output.matrix[(i*2),k] <- "" #lower std row
        } else if (m[i,j] == -Inf) {
          output.matrix[(i*2)-1,k] <- neginfstring #upper row
          output.matrix[(i*2),k] <- "" #lower std row
        } else {
          if (strong.signif == TRUE && stars==TRUE) {
            if (m[i,j+2] <= 0.001) {
              p <- paste(star.prefix, "***", star.suffix, sep="")
            } else if (m[i,j+2] <= 0.01) {
              p <- paste(star.prefix, "**", star.suffix, sep="")
            } else if (m[i,j+2] <= 0.05) {
              p <- paste(star.prefix, "*", star.suffix, sep="")
            } else if (m[i,j+2] <= 0.1) {
              p <- paste(star.prefix, symbol, star.suffix, sep="")
            } else {
              p <- ""
            }
          } else if (stars==TRUE) {
            if (m[i,j+2] <= 0.01) {
              p <- paste(star.prefix, "***", star.suffix, sep="")
            } else if (m[i,j+2] <= 0.05) {
              p <- paste(star.prefix, "**", star.suffix, sep="")
            } else if (m[i,j+2] <= 0.1) {
              p <- paste(star.prefix, "*", star.suffix, sep="")
            } else {
              p <- ""
            }
          } else {
              p <- ""
          }
          if (dcolumn == TRUE) {
            dollar <- ""
          } else {
            dollar <- "$"
          }
          output.matrix[(i*2)-1,k] <- paste(dollar, coeftostring(m[i,j], 
              leading.zero, digits=digits), p, dollar, sep="")
          output.matrix[(i*2),k] <- paste(dollar, "(", 
              coeftostring(m[i,j+1], leading.zero, digits=digits), ")", 
              dollar, sep="")
        }
        k <- k+1
        j <- j+3
      }
    }
  }
  
  return(output.matrix)
}


# Format a column (given as vector) of the output matrix nicely by adding spaces
format.column <- function(x, single.row=FALSE, digits=2) {
  
  #max length before first dot and max length of parentheses
  dots <- gregexpr("\\.", x)
  parentheses <- regexpr("\\(.+\\)", x)
  first.length <- 0
  paren.length <- 0
  for (i in 1:length(x)) {
    first.dot <- dots[[i]][1]
    paren <- attributes(parentheses)$match.length[i]
    if (x[i] == "-Inf") {
      first.dot <- nchar(x[i]) - digits
    } else if (first.dot == -1) {
      temp <- nchar(x[i]) + 1
      if (temp > first.length) {
        first.length <- temp
      }
    } else if (first.dot > first.length) {
      first.length <- first.dot
    }
    if (paren > paren.length) {
      paren.length <- paren
    }
  }
  
  for (i in 1:length(x)) {
    
    #fill with spaces at the beginning
    first.dot <- dots[[i]][1]
    if (x[i] == "-Inf") {
      first.dot <- nchar(x[i]) - digits
    } else if (first.dot == -1) {
      first.dot <- nchar(x[i]) + 1
    }
    if (nchar(x[i]) == 0) {
      difference <- 0
    } else {
      difference <- first.length - first.dot
    }
    spaces <- paste(rep(" ", difference), collapse="")
    x[i] <- paste(spaces, x[i], sep="")
    
    #adjust indentation for SEs
    if (single.row==TRUE) {
      paren <- attributes(parentheses)$match.length[i]
      if (paren < 0) {
        paren <- 0
      }
      difference <- paren.length - paren + 1 #+1 because strsplit takes one away
      spaces <- paste(rep(" ", difference), collapse="")
      components <- strsplit(x[i], " \\(")[[1]]
      if (length(components) == 2) {
        x[i] <- paste(components[1], spaces, "(", components[2], sep="")
      }
    }
  }
  
  #fill with spaces at the end to make them all equally long
  max.x <- max(nchar(x))
  for (i in 1:length(x)) {
    difference <- max.x - nchar(x[i])
    spaces <- paste(rep(" ", difference), collapse="")
    x[i] <- paste(x[i], spaces, sep="")
  }

  return(x)
}


# fill a column/vector with spaces at the end
fill.spaces <- function(x) {
  nc <- nchar(x)
  width <- max(nc)
  for (i in 1:length(x)) {
    spaces <- paste(rep(" ", width - nc[i]), collapse="")
    x[i] <- paste(x[i], spaces, sep="")
  }
  return(x)
}


# Return the goodness-of-fit matrix (i.e., the lower block of the final matrix)
gofmatrix <- function(gofs, decimal.matrix, dcolumn=TRUE, leading.zero, 
    digits) {
  if (dcolumn == TRUE) {
    dollar <- ""
  } else {
    dollar <- "$"
  }
  gof.matrix <- matrix(nrow=nrow(gofs), ncol=ncol(gofs)+1) #incl. labels
  for (i in 1:length(gofs[,1])) {
    gof.matrix[i,1] <- rownames(gofs)[i]
    for (j in 1:length(gofs[1,])) {
      strg <- coeftostring(gofs[i,j], leading.zero, digits=decimal.matrix[i,j])
      gof.matrix[i,j+1] <- paste(dollar, strg, dollar, sep="")
    }
  }
  return(gof.matrix)
}


