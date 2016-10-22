# The texreg package was written by Philip Leifeld.
# Please use the forum at http://r-forge.r-project.org/projects/texreg/ 
# for bug reports, help or feature requests.


# display version number and date when the package is loaded
.onAttach <- function(libname, pkgname) {
  desc  <- packageDescription(pkgname, libname)
  packageStartupMessage(
      'Version:  ', desc$Version, '\n', 
      'Date:     ', desc$Date, '\n',
      'Author:   ', 'Philip Leifeld (University of Glasgow)', '\n\n', 
      'Please cite the JSS article in your publications ', 
      '-- see citation("texreg").'
  )
}


# function which reformats a coefficient with two decimal places
coeftostring <- function(x, lead.zero = FALSE, digits = 2) {
  if (!is.finite(x)) {
    return("")
  }
  if (digits < 0) {
    stop("The number of digits must be 0 or higher.")
  }
  
  y <- format(round(x, digits), nsmall = digits, scientific = FALSE)
  
  if (lead.zero == FALSE && (grepl("^0", y) == TRUE ||  # leading zero
      grepl("^-0", y) == TRUE)) {
    y <- gsub("0\\.", "\\.", y)
  }
  if (x < 0 && grepl("^-", y) == FALSE) {  # very small negative numbers
    y <- paste0("-", y)
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
  orig.width <- length(m[1, ])                     #number of columns in m
  q <- matrix(nrow = 0, ncol = orig.width)         #new matrix with same width
  for (i in 1:num.unique) {                        #go through unique row names
    rows <- matrix(NA, nrow = 0, ncol = orig.width)#create matrix where re-
                                                   #arranged rows will be stored
    for (j in 1:orig.width) {                      #go through columns in m
      current.name <- unique.names[i]              #save row name
      nonNa <- m[rownames(m) == current.name, j]   #create a vector of values
                                                   #with same rowname in the col
      nonNa <- nonNa[!is.na(nonNa)]                #retain only non-NA values
      for (k in 1:length(nonNa)) {                 #go through non-NA values
        if (k > dim(rows)[1]) {                    #add an NA-only row in which
          rows <- rbind(rows, rep(NA, orig.width)) #the values are stored
          rownames(rows)[k] <- unique.names[i]     #also add the row name
        }
        rows[k, j] <- nonNa[k]                     #actually store the value
      }
    }
    q <- rbind(q, rows)                            #add the new row(s) to q
  }
  return(q)
}


# function which wraps models in a list and extracts texreg objects from them
get.data <- function(l, ...) {

  # if a single model is handed over, put model inside a list
  if (!"list" %in% class(l)[1]) {
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
    if (!is.null(gn) && length(gn) > 0) {
      for (j in 1:length(gn)) {
        if (!gn[j] %in% gof.names) {
          gof.names <- append(gof.names, gn[j])
        }
      }
    }
  }
  return(gof.names)
}


# function which replaces coefs, SEs and p values by custom values if provided
override <- function(models, override.coef, override.se, override.pval, 
    override.ci.low, override.ci.up) {
  
  if (class(override.se) == "list" || length(override.se) > 1 || 
      override.se[1] != 0) {
    if (length(override.pval) == 1 && class(override.pval) != "list" && 
        override.pval[1] == 0) {
      warning(paste("Standard errors were provided using 'override.se',", 
          "but p-values were not replaced!"))
    }
  }
  if (class(override.pval) == "list" || length(override.pval) > 1 || 
      override.pval[1] != 0) {
    if (length(override.se) == 1 && class(override.se) != "list" && 
        override.se[1] == 0) {
      warning(paste("P-values were provided using 'override.pval',", 
          "but standard errors were not replaced!"))
    }
  }
  
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
          "coefficients in model ", i, ". Using default SEs.", sep = ""))
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
          "coefficients in model ", i, ". Using default p values.", sep = ""))
      pval <- models[[i]]@pvalues
    } else if (class(override.pval[[i]]) != "numeric") {
      warning(paste("p values provided for model", i, 
          "are not numeric. Using default p values."))
      pval <- models[[i]]@pvalues
    } else {
      pval <- override.pval[[i]]
    }
    models[[i]]@pvalues <- pval
    
    # lower bound of confidence intervals
    if (is.null(override.ci.low)) {
      # do nothing
    } else if (class(override.ci.low) != "list" && length(override.ci.low) 
        == 1 && override.ci.low == 0) {
      ci.low <- models[[i]]@ci.low
    } else if (class(override.ci.low) == "numeric" && length(models) == 1 && 
        length(override.ci.low) == length(models[[i]]@coef)) {
      ci.low <- override.ci.low
    } else if (class(override.ci.low) != "list") {
      warning("CIs must be provided as a list. Using default CIs if available.")
      ci.low <- models[[i]]@ci.low
    } else if (length(override.ci.low) != length(models)) {
      warning(paste("Number of lower CIs provided does not match number of", 
          "models. Using default CIs if available."))
      ci.low <- models[[i]]@ci.low
    } else if (length(models[[i]]@coef) != length(override.ci.low[[i]])) {
      # previous line: comparison with coef because CIs can be empty
      warning(paste0("Number of lower CIs provided does not match number of ", 
          "coefficients in model ", i, ". Using default CIs if available."))
      ci.low <- models[[i]]@ci.low
    } else if (class(override.ci.low[[i]]) != "numeric") {
      warning(paste("Lower CIs provided for model", i, 
          "are not numeric. Using default lower CIs."))
      ci.low <- models[[i]]@ci.low
    } else {
      ci.low <- override.ci.low[[i]]
    }
    models[[i]]@ci.low <- ci.low
    
    # upper bound of confidence intervals
    if (is.null(override.ci.up)) {
      # do nothing
    } else if (class(override.ci.up) != "list" && length(override.ci.up) 
        == 1 && override.ci.up == 0) {
      ci.up <- models[[i]]@ci.up
    } else if (class(override.ci.up) == "numeric" && length(models) == 1 && 
        length(override.ci.up) == length(models[[i]]@coef)) {
      ci.up <- override.ci.up
    } else if (class(override.ci.up) != "list") {
      warning("CIs must be provided as a list. Using default CIs if available.")
      ci.up <- models[[i]]@ci.up
    } else if (length(override.ci.up) != length(models)) {
      warning(paste("Number of lower CIs provided does not match number of", 
          "models. Using default CIs if available."))
      ci.up <- models[[i]]@ci.up
    } else if (length(models[[i]]@coef) != length(override.ci.up[[i]])) {
      # previous line: comparison with coef because CIs can be empty
      warning(paste0("Number of lower CIs provided does not match number of ", 
          "coefficients in model ", i, ". Using default CIs if available."))
      ci.up <- models[[i]]@ci.up
    } else if (class(override.ci.up[[i]]) != "numeric") {
      warning(paste("Lower CIs provided for model", i, 
          "are not numeric. Using default lower CIs."))
      ci.up <- models[[i]]@ci.up
    } else {
      ci.up <- override.ci.up[[i]]
    }
    models[[i]]@ci.up <- ci.up
    
    if (length(models[[i]]@ci.low) > 0 && length(models[[i]]@ci.up) > 0) {
      models[[i]]@se <- numeric()
      models[[i]]@pvalues <- numeric()
    }
  }
  
  return(models)
}


# function which converts GOF/coef LaTeX code to HTML oder text/screen code
tex.replace <- function(models, type = "html", style = "") {
  for (i in 1:length(models)) {
    # GOF replacement
    if (type == "html") {
      r <- paste0("<sup", style, ">2</sup>")
    } else if (type == "screen") {
      r <- "^2"
    }
    models[[i]]@gof.names <- gsub("\\$\\^2\\$", r, models[[i]]@gof.names)
    models[[i]]@gof.names <- gsub("\\\\ ", " ", models[[i]]@gof.names)
    models[[i]]@gof.names <- gsub("\\ ", " ", models[[i]]@gof.names)
    
    # extract.sarlm coefficient name replacement
    models[[i]]@coef.names <- gsub("\\$\\\\rho\\$", "rho", models[[i]]@coef.names)
    models[[i]]@coef.names <- gsub("\\$\\\\lambda\\$", "lambda", models[[i]]@coef.names)
    
    # extract.gamlss coefficient name replacement
    models[[i]]@coef.names <- gsub("\\$\\\\mu\\$", "mu", models[[i]]@coef.names)
    models[[i]]@coef.names <- gsub("\\$\\\\nu\\$", "nu", models[[i]]@coef.names)
    models[[i]]@coef.names <- gsub("\\$\\\\tau\\$", "tau", 
        models[[i]]@coef.names)
    models[[i]]@coef.names <- gsub("\\$\\\\sigma\\$", "sigma", 
        models[[i]]@coef.names)
  }
  return(models)
}


# function which replaces special characters in row names by LaTeX equivalents
replaceSymbols <- function(m) {
  rn <- rownames(m)
  for (i in 1:length(rn)) {
    if (!grepl("\\$", rn[i])) {
      rn[i] <- gsub("_", "\\\\_", rn[i])
      rn[i] <- gsub("<", "\\$<\\$", rn[i])
      rn[i] <- gsub(">", "\\$>\\$", rn[i])
      rn[i] <- gsub("%", "\\\\%", rn[i])
    }
  }
  rownames(m) <- rn
  return(m)
}


# correct duplicate coefficient names (add " (1)", " (2)" etc.)
correctDuplicateCoefNames <- function(models) {
  for (i in 1:length(models)) {
    for (j in 1:length(models[[i]]@coef.names)) {
      if (models[[i]]@coef.names[j] %in% models[[i]]@coef.names[-j]) {
        indices <- j
        for (k in 1:length(models[[i]]@coef.names)) {
          if (models[[i]]@coef.names[j] == models[[i]]@coef.names[k] && j != k) {
            indices <- c(indices, k)
          }
        }
        count <- 1
        for (k in indices) {
          models[[i]]@coef.names[k] <- paste0(models[[i]]@coef.names[k], " (", 
              count, ")")
          count <- count + 1
        }
      }
    }
  }
  return(models)
}


# put models and GOFs into a common matrix
aggregate.matrix <- function(models, gof.names, custom.gof.names, digits, 
    returnobject = "m") {

  # aggregate GOF statistics in a matrix and create list of coef blocks
  gofs <- matrix(nrow = length(gof.names), ncol = length(models))
  row.names(gofs) <- gof.names
  coefs <- list()
  decimal.matrix <- matrix(nrow = length(gof.names), ncol = length(models))
  for (i in 1:length(models)) {
    cf <- models[[i]]@coef
    se <- models[[i]]@se
    pv <- models[[i]]@pvalues
    cil <- models[[i]]@ci.low
    ciu <- models[[i]]@ci.up
    if (length(se) == 0 && length(ciu) > 0) {
      coef <- cbind(cf, cil, ciu)
    } else {
      if (length(se) > 0 && length(pv) > 0) {
        coef <- cbind(cf, se, pv)
      } else if (length(se) > 0 && length(pv) == 0) {
        #p-values not provided -> use p-values of 0.99
        coef <- cbind(cf, se, rep(0.99, length(cf)))
      } else if (length(se) == 0 && length(pv) > 0) {
        coef <- cbind(cf, rep(NA, length(cf)), pv)
      } else {
        # not even SEs provided
        coef <- cbind(cf, rep(NA, length(cf)), rep(0.99, length(cf)))
      }
    }
    rownames(coef) <- models[[i]]@coef.names
    coefs[[i]] <- coef
    if (length(models[[i]]@gof) > 0) {
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
        gofs[row, col] <- val
        decimal.matrix[row, col] <- dec
      }
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
      m <- merge(m, coefs[[i]], by = 0, all = TRUE)
      rownames(m) <- m[, 1]
      m <- m[, colnames(m) != "Row.names"]
      colnames(m) <- NULL
    }
  }
  colnames(m) <- rep(colnames(coefs[[1]]), length(coefs))
  
  # reorder merged coefficient table
  m.temp <- matrix(nrow = nrow(m), ncol = ncol(m))
  for (i in 1:nrow(m)) {
    new.row <- which(coef.order == rownames(m)[i])
    for (j in 1:length(m[i,])) {
      m.temp[new.row, j] <- m[i, j]
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
    } else {
      custom.gof.names[is.na(custom.gof.names)] <- 
          rownames(gofs)[is.na(custom.gof.names)]
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
      custom.names[is.na(custom.names)] <- rownames(m)[is.na(custom.names)]
      rownames(m) <- custom.names
    }
  } else if (!is.na(custom.names) & class(custom.names) != "character") {
    stop("Custom coefficient names must be provided as a vector of strings.")
  } else if (length(custom.names) == 1 & class(custom.names) == "character"
      & is.na(custom.names)) {
    rownames(m) <- custom.names
  } else if (length(custom.names) == length(rownames(m))) {
    rownames(m) <- custom.names
  }
  return(m)
}


# remove coefficient rows that match the omit.coef regular expression
omitcoef <- function(m, omit.coef) {
  if (!is.null(omit.coef)) {
    if (!is.character(omit.coef) || is.na(omit.coef)) {
      stop("omit.coef must be a character string!")
    }
    remove.rows <- grep(omit.coef, rownames(m))
    if (length(remove.rows) == 0) {
      return(m)
    } else if (length(remove.rows) == nrow(m)) {
      stop("You were trying to remove all coefficients using omit.coef.")
    } else {
      m <- m[-remove.rows, ]
    }
  }
  return(m)
}


# decide if default or custom model names should be used and return them
modelnames <- function(model.list, tr.objects, model.names) {
  
  if (class(model.list)[1] != "list") {
    model.list <- list(model.list)
  }
  
  mnames <- names(model.list)
  if (is.null(mnames)) {
    mnames <- character(length(model.list))
    if (length(mnames) != length(tr.objects)) {
      mnames <- character(length(tr.objects))
    }
  }
  
  for (i in 1:length(tr.objects)) {
    nam <- tr.objects[[i]]@model.name
    if (length(nam) == 1) {
      model.names[i] <- nam
    }
  }
  
  if (is.null(model.names)) {
    model.names <- rep(NA, length(mnames))
  } else if (class(model.names) != "character") {
    stop("Model names must be specified as a vector of strings.")
  } else if (length(model.names) != length(tr.objects)) {
    stop(paste("There are", length(tr.objects), "models, but you provided", 
        length(model.names), "name(s) for them."))
  }
  
  for (i in 1:length(model.names)) {
    if (!is.na(model.names[i])) {
      mnames[i] <- model.names[i]
    } else if (mnames[i] == "") {
      mnames[i] <- paste("Model", i)
    } else {
      # leave mnames[i] as is
    }
  }
  
  return(mnames)
}


# check if stars argument is OK
check.stars <- function(stars) {
  if (is.null(stars)) {
    s <- numeric(0)
  } else if (class(stars) != "numeric") {
    stop("The stars argument must be a numeric vector.")
  } else if (any(is.na(stars))) {
    stop("NA value are not allowed in the stars argument.")
  } else {
    s <- stars
  }
  return(s)
}


# create stars string
stars.string <- function(pval, stars, star.char, star.prefix, star.suffix, 
    symbol) {
  st <- sort(stars)
  if (length(unique(st)) != length(st)) {
    stop("Duplicate elements are not allowed in the stars argument.")
  }
  if (is.na(pval)) {
    pval <- 1.0
  }
  if (length(st) > 4) {
    stop("A maximum of four values is allowed in the stars argument.")
  } else if (length(st) > 2 && pval < st[1]) {  # three stars
    p <- paste0(star.prefix, star.char, star.char, star.char, star.suffix)
  } else if (  # two stars
      (length(st) > 2 && pval < st[2]) || 
      (length(st) == 2 && pval < st[1]) ) {
    p <- paste0(star.prefix, star.char, star.char, star.suffix)
  } else if (  # one star
      (length(st) > 2 && pval < st[3]) || 
      (length(st) == 2 && pval < st[2]) || 
      (length(st) == 1 && pval < st) ) {
    p <- paste0(star.prefix, star.char, star.suffix)
  } else if (length(st) == 4 && pval < st[4]) {  # symbol
    p <- paste0(star.prefix, symbol, star.suffix)
  } else {  # not significant
    p <- ""
  }
  return(p) 
}


# return the output matrix with coefficients, SEs and significance stars
outputmatrix <- function(m, single.row, neginfstring, posinfstring, 
    leading.zero, digits, se.prefix, se.suffix, star.prefix, star.suffix, 
    star.char = "*", stars, dcolumn = TRUE, symbol, bold, bold.prefix, 
    bold.suffix, ci = rep(FALSE, length(m) / 3), semicolon = "; ", 
    ci.test = 0) {
  
  # write coefficient rows
  if (single.row == TRUE) {
    output.matrix <- matrix(ncol = (length(m) / 3) + 1, nrow = length(m[, 1]))
    
    # row labels
    for (i in 1:length(rownames(m))) {
      output.matrix[i, 1] <- rownames(m)[i]
    }
    
    # coefficients and standard errors
    for (i in 1:length(m[, 1])) { #go through rows
      j <- 1 #column in the original, merged coef table
      k <- 2 #second column of output.matrix, i.e., coefficients
      while (j <= length(m)) {
        if (is.na(m[i, j])) {
          output.matrix[i, k] <- ""
        } else if (m[i, j] == -Inf) {
          output.matrix[i, k] <- neginfstring
        } else if (m[i, j] == Inf) {
          output.matrix[i, k] <- posinfstring
        } else {
          
          # in case of CIs, replace brackets by square brackets
          se.prefix.current <- se.prefix
          se.suffix.current <- se.suffix
          if (ci[k - 1] == TRUE) {
            se.prefix.current <- gsub("\\(", "[", se.prefix.current)
            se.suffix.current <- gsub("\\)", "]", se.suffix.current)
          }
          if (is.na(m[i, j + 1])) {
            se.prefix.current <- ""
            se.suffix.current <- ""
          }
          if (ci[k - 1] == FALSE) {
            std <- paste(se.prefix.current, coeftostring(m[i, j + 1], 
                leading.zero, digits = digits), se.suffix.current, sep = "")
          } else {
            std <- paste(se.prefix.current, coeftostring(m[i, j + 1], 
                leading.zero, digits = digits), semicolon, 
                coeftostring(m[i, j + 2], leading.zero, digits = digits), 
                se.suffix.current, sep = "")
          }
          
          if (ci[k - 1] == FALSE) {
            p <- stars.string(m[i, j + 2], stars, star.char, star.prefix, 
              star.suffix, symbol)
          } else { # significance from confidence interval
            if (is.numeric(ci.test) && !is.na(ci.test) && 
                (m[i, j + 1] > ci.test || m[i, j + 2] < ci.test)) {
              p <- paste0(star.prefix, star.char, star.suffix)
            } else {
              p <- ""
            }
          }
          
          if (dcolumn == TRUE) {
            dollar <- ""
          } else {
            dollar <- "$"
          }
          if (is.na(m[i, j + 2])) {
            m[i, j + 2] <- 1.0
          }
          if (ci[k - 1] == FALSE && m[i, j + 2] < bold) { # significant p-value
            bold.pref <- bold.prefix
            bold.suff <- bold.suffix
          } else if (ci[k - 1] == TRUE && bold > 0 &&  # significant CI
              (m[i, j + 1] > 0 || m[i, j + 2] < 0)) {
            bold.pref <- bold.prefix
            bold.suff <- bold.suffix
          } else {
            bold.pref <- ""
            bold.suff <- ""
          }
          entry <- paste(dollar, bold.pref, coeftostring(m[i, j], leading.zero, 
              digits = digits), bold.suff, std, p, dollar, sep = "")
          output.matrix[i, k] <- entry
          
        }
        k <- k + 1
        j <- j + 3
      }
    }
  } else {
    output.matrix <- matrix(ncol = (length(m) / 3) + 1, 
        nrow = 2 * length(m[, 1]))
    
    # row labels
    for (i in 1:length(rownames(m))) {
      output.matrix[(i * 2) - 1, 1] <- rownames(m)[i]
      output.matrix[(i * 2), 1] <- ""
    }
    
    # coefficients and standard deviations
    for (i in 1:length(m[, 1])) {  # i = row
      j <- 1  # j = column within model (from 1 to 3)
      k <- 2  # k = column in output matrix (= model number + 1)
      while (j <= length(m)) {
        if (is.na(m[i, j]) || is.nan(m[i, j])) {
          output.matrix[(i * 2) - 1, k] <- ""  #upper coefficient row
          output.matrix[(i * 2), k] <- ""  #lower se row
        } else if (m[i, j] == -Inf) {
          output.matrix[(i * 2) - 1, k] <- neginfstring  #upper row
          output.matrix[(i * 2), k] <- ""  #lower se row
        } else if (m[i, j] == Inf) {
          output.matrix[(i * 2) - 1, k] <- posinfstring  #upper row
          output.matrix[(i * 2), k] <- ""  #lower se row
        } else {
          
          # in case of CIs, replace brackets by square brackets
          se.prefix.current <- "("
          se.suffix.current <- ")"
          if (ci[k - 1] == TRUE) {
            se.prefix.current <- "["
            se.suffix.current <- "]"
          }
          if (is.na(m[i, j + 1])) {
            se.prefix.current <- ""
            se.suffix.current <- ""
          }
          if (ci[k - 1] == FALSE) {
            p <- stars.string(m[i, j + 2], stars, star.char, star.prefix, 
              star.suffix, symbol)
          } else { # significance from confidence interval
            if (is.numeric(ci.test) && !is.na(ci.test) && 
                (m[i, j + 1] > ci.test || m[i, j + 2] < ci.test)) {
              p <- paste0(star.prefix, star.char, star.suffix)
            } else {
              p <- ""
            }
          }
          
          if (dcolumn == TRUE) {
            dollar <- ""
          } else {
            dollar <- "$"
          }
          if (is.na(m[i, j + 2])) {
            m[i, j + 2] <- 1.0
          }
          if (ci[k - 1] == FALSE && m[i, j + 2] < bold) { # significant p-value
            bold.pref <- bold.prefix
            bold.suff <- bold.suffix
          } else if (ci[k - 1] == TRUE && bold > 0 &&  # significant CI
              (m[i, j + 1] > 0 || m[i, j + 2] < 0)) {
            bold.pref <- bold.prefix
            bold.suff <- bold.suffix
          } else {
            bold.pref <- ""
            bold.suff <- ""
          }
          output.matrix[(i * 2) - 1, k] <- paste(dollar, bold.pref, 
              coeftostring(m[i, j], leading.zero, digits = digits), bold.suff, 
              p, dollar, sep = "")
          if (ci[k - 1] == FALSE) {
            output.matrix[(i * 2), k] <- paste(dollar, se.prefix.current, 
                coeftostring(m[i, j + 1], leading.zero, digits = digits), 
                se.suffix.current, dollar, sep = "")
          } else {
            output.matrix[(i * 2), k] <- paste(dollar, se.prefix.current, 
                coeftostring(m[i, j + 1], leading.zero, digits = digits), 
                semicolon, coeftostring(m[i, j + 2], leading.zero, 
                digits = digits), se.suffix.current, dollar, sep = "")
          }
        }
        k <- k + 1
        j <- j + 3
      }
    }
    
    # check if SEs are all missing and delete even rows if necessary
    se.missing <- numeric()
    for (i in seq(2, nrow(output.matrix), 2)) {
      if (all(sapply(output.matrix[i, ], function(x) x == ""))) {
        se.missing <- c(se.missing, i)
      }
    }
    if (length(se.missing) == nrow(output.matrix) / 2) {
      output.matrix <- output.matrix[-se.missing, ]
    }
  }
  
  return(output.matrix)
}


# Format a column (given as vector) of the output matrix nicely by adding spaces
format.column <- function(x, single.row = FALSE, digits = 2) {
  
  #max length before first dot and max length of parentheses
  dots <- gregexpr("\\.", x)
  parentheses <- regexpr("\\(.+\\)", x)
  first.length <- 0
  paren.length <- 0
  for (i in 1:length(x)) {
    first.dot <- dots[[i]][1]
    paren <- attributes(parentheses)$match.length[i]
    if (x[i] %in% c("-Inf", "Inf")) {
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
    if (x[i] %in% c("-Inf", "Inf")) {
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
    if (single.row == TRUE) {
      paren <- attributes(parentheses)$match.length[i]
      if (paren < 0) {
        paren <- 0
      }
      difference <- paren.length - paren + 1  #+1 because strsplit takes 1 away
      spaces <- paste(rep(" ", difference), collapse = "")
      components <- strsplit(x[i], " \\(")[[1]]
      if (length(components) == 2) {
        x[i] <- paste(components[1], spaces, "(", components[2], sep = "")
      }
    }
  }
  
  #make all CIs have equal length
  ci.lower.length <- 0
  ci.upper.length <- 0
  for (i in 1:length(x)) {
    if (grepl("\\[.+\\]", x[i])) {
      first <- sub(".*\\[(.+?); (.+?)\\].*", "\\1", x[i])
      first <- nchar(first)
      if (first > ci.lower.length) {
        ci.lower.length <- first
      }
      last <- sub(".*\\[(.+?); (.+?)\\].*", "\\2", x[i])
      last <- nchar(last)
      if (last > ci.upper.length) {
        ci.upper.length <- last
      }
    }
  }
  for (i in 1:length(x)) {
    if (grepl("\\[.+\\]", x[i])) {
      whitespace1 <- sub("(.*?)\\[(.+?); (.+?)\\](.*?)$", "\\1", x[i])
      whitespace1 <- sub("\\s+$", "", whitespace1)
      if (nchar(whitespace1) > 0) {
        whitespace1 <- paste0(whitespace1, " ")
      }
      whitespace2 <- sub("(.*?)\\[(.+?); (.+?)\\](.*?)$", "\\4", x[i])
      first <- sub("(.*?)\\[(.+?); (.+?)\\](.*?)$", "\\2", x[i])
      difference <- ci.lower.length - nchar(first)
      zeros <- paste(rep(" ", difference), collapse = "")
      first <- paste0(zeros, first)
      last <- sub("(.*?)\\[(.+?); (.+?)\\](.*?)$", "\\3", x[i])
      difference <- ci.upper.length - nchar(last)
      zeros <- paste(rep(" ", difference), collapse = "")
      last <- paste0(zeros, last)
      x[i] <- paste0(whitespace1, "[", first, "; ", last, "]", whitespace2)
    }
  }
  
  #fill with spaces at the end to make them all equally long
  max.x <- max(nchar(x))
  for (i in 1:length(x)) {
    difference <- max.x - nchar(x[i])
    spaces <- paste(rep(" ", difference), collapse = "")
    x[i] <- paste(x[i], spaces, sep = "")
  }

  return(x)
}


# fill a column/vector with spaces at the end
fill.spaces <- function(x) {
  nc <- nchar(x)
  width <- max(nc)
  for (i in 1:length(x)) {
    spaces <- paste(rep(" ", width - nc[i]), collapse = "")
    x[i] <- paste(x[i], spaces, sep = "")
  }
  return(x)
}


# Return the goodness-of-fit matrix (i.e., the lower block of the final matrix)
gofmatrix <- function(gofs, decimal.matrix, dcolumn = TRUE, leading.zero, 
    digits) {
  if (dcolumn == TRUE) {
    dollar <- ""
  } else {
    dollar <- "$"
  }
  gof.matrix <- matrix(nrow = nrow(gofs), ncol = ncol(gofs) + 1)  #incl. labels
  if (length(gof.matrix) > 0) {
    for (i in 1:length(gofs[, 1])) {
      gof.matrix[i, 1] <- rownames(gofs)[i]
      for (j in 1:length(gofs[1, ])) {
        strg <- coeftostring(gofs[i, j], leading.zero, 
            digits = decimal.matrix[i, j])
        gof.matrix[i, j + 1] <- paste0(dollar, strg, dollar)
      }
    }
  }
  return(gof.matrix)
}


# reorder a matrix according to a vector of new positions
reorder <- function(mat, new.order) {
  if (is.null(new.order)) {
    return(mat)
  } else if (nrow(mat) != length(new.order)) {
    stop(paste("Error when reordering matrix: there are", nrow(mat), 
        "rows, but you provided", length(new.order), "numbers."))
  } else if (class(new.order) == "list") {
    stop("Arguments reorder.coef and reorder.gof must be provided as a vector.")
  } else if (any(is.na(new.order))) {
    stop("reorder.coef and reorder.gof arguments must not contain NA values.")
  } else if (length(new.order) != length(unique(new.order))) {
    stop(paste("There are two identical values in the reorder.coef or", 
        "reorder.gof argument. Ties are not allowed."))
  } else if (max(new.order) != nrow(mat)) {
    stop(paste("Table cannot be reordered because you provided a number that",
        "exceeds the number of rows of the relevant part of the table."))
  }
  new.sorted <- sort(new.order)
  for (i in 2:length(new.sorted)) {
    if (new.sorted[i] - 1 != new.sorted[i - 1]) {
      stop(paste("Table cannot be reordered because there are non-adjacent", 
          "values in the reorder.coef or reorder.gof vector you provided."))
    }
  }
  new.mat <- mat[new.order, ]
  return(new.mat)
}


# compute column width left and right of the decimal separator
compute.width <- function(v, left = TRUE, single.row = FALSE, bracket = ")") {
  if (single.row == FALSE) {
    v[which(!grepl("\\.", v))] <- paste0(v[which(!grepl("\\.", v))], ".")
    ssp <- strsplit(v, "\\.")
    left.side <- character()
    right.side <- character()
    for (i in 1:length(ssp)) {
      if (length(ssp[[i]]) == 1) {
        ssp[[i]][2] <- ""
      } else if (length(ssp[[i]]) == 3) {
        ssp[[i]] <- c(ssp[[i]][1], paste0(ssp[[i]][2], ".", ssp[[i]][3]))
      }
      left.side[i] <- ssp[[i]][1]
      right.side[i] <- ssp[[i]][2]
    }
  } else {
    ssp <- strsplit(v, paste0("\\", bracket))
    left.side <- character()
    right.side <- character()
    for (i in 1:length(ssp)) {
      if (length(ssp[[i]]) == 0) {
        # do nothing because empty cell
      } else {
        left.side <- append(left.side, ssp[[i]][1])
        right.side <- append(right.side, ssp[[i]][2])
      }
    }
  }
  if (left == TRUE) {
    left.side <- sub("\\\\; ", "", left.side)
    v.length <- max(nchar(left.side), na.rm = TRUE)
  } else {
    right.side <- sub("\\^\\{", "", right.side)
    right.side <- sub("\\}", "", right.side)
    v.length <- max(nchar(right.side), na.rm = TRUE)
  }
  return(v.length)
}


# convert SEs and p values to confidence intervals
ciforce <- function(models, ci.force = rep(FALSE, length(models)), 
    ci.level = 0.95) {
  if (class(ci.force) == "logical" && length(ci.force) == 1) {
    ci.force <- rep(ci.force, length(models))
  }
  if (class(ci.force) != "logical") {
    stop("The 'ci.force' argument must be a vector of logical values.")
  }
  if (length(ci.force) != length(models)) {
    stop(paste("There are", length(models), "models and", length(ci.force), 
        "ci.force values."))
  }
  for (i in 1:length(models)) {
    if (ci.force[i] == TRUE && length(models[[i]]@se) > 0) {
      z <- qnorm(1 - ((1 - ci.level) / 2))
      upper <- models[[i]]@coef + (z * models[[i]]@se)
      lower <- models[[i]]@coef - (z * models[[i]]@se)
      models[[i]]@ci.low <- lower
      models[[i]]@ci.up <- upper
      models[[i]]@se <- numeric(0)
      models[[i]]@pvalues <- numeric(0)
    }
  }
  return(models)
}

# function which adds groups to an output matrix
grouping <- function(output.matrix, groups, indentation = "    ", 
    single.row = FALSE, prefix = "", suffix = "") {
  if (!is.null(groups)) {
    if (class(groups) != "list") {
      stop("Groups must be specified as a list of numeric vectors.")
    }
    for (i in 1:length(groups)) {
      if (length(groups[[i]]) == 0) {
        stop("Empty groups are not allowed.")
      }
      if (!is.numeric(groups[[i]])) {
        stop("Groups must be specified as a list of numeric vectors.")
      }
      groups[[i]] <- sort(unique(groups[[i]]))
      if (groups[[i]][length(groups[[i]])] - length(groups[[i]]) + 1 != 
          groups[[i]][1]) {
        stop("The group indices must be consecutive within each group.")
      }
    }
    for (i in 1:length(groups)) {
      for (j in 1:length(groups)) {
        if (i != j && length(intersect(groups[[i]], groups[[j]])) > 0) {
          stop("Overlapping groups are not allowed. Change 'groups' argument!")
        }
        if (i < j && groups[[j]][1] < groups[[i]][1]) {
          stop("Groups must be specified in the correct order.")
        }
      }
    }
    for (i in 1:length(groups)) {
      if (single.row == FALSE) {
        groups[[i]] <- (groups[[i]] * 2) - 1
      }
    }
    if (groups[[length(groups)]][length(groups[[length(groups)]])] > 
        nrow(output.matrix)) {
      stop("'groups' argument contains indices outside the table dimensions")
    }
    for (i in length(groups):1) {
      label <- paste0(prefix, names(groups)[[i]], suffix)
      for (j in nrow(output.matrix):1) {
        if (j %in% groups[[i]]) {
          output.matrix[j, 1] <- paste0(indentation, output.matrix[j, 1])
        }
      }
      groupindex <- groups[[i]][1]
      lastingroup <- groups[[i]][length(groups[[i]])] + (1 - single.row)
      if (groupindex == 1) {
        prevmat <- NULL
      } else {
        prevmat <- output.matrix[1:(groupindex - 1), ]
      }
      currentmat <- output.matrix[groupindex:lastingroup, ]
      if (lastingroup > nrow(output.matrix) - 2) {
        nextmat <- NULL
      } else {
        nextmat <- output.matrix[(lastingroup + 1):nrow(output.matrix), ]
      }
      newrow <- matrix(rep("", ncol(output.matrix)), nrow = 1)
      if (single.row == FALSE) {
        newrow <- rbind(newrow, newrow)
      }
      newrow[1, 1] <- label
      output.matrix <- rbind(prevmat, newrow, currentmat, nextmat)
    }
  }
  return(output.matrix)
}

# add custom columns to output matrix
customcolumns <- function(output.matrix, custom.columns, custom.col.pos, 
    single.row = FALSE, numcoef, groups, modelnames = TRUE) {
  
  # check validity of arguments
  if (is.null(custom.columns)) {
    return(output.matrix)
  }
  if (!class(custom.columns) == "list") {
    if (length(custom.columns) != numcoef) {
      stop(paste("Custom column does not match table dimensions.", numcoef, 
          "elements expected."))
    }
    custom.columns <- list(custom.columns)
  }
  if (is.null(custom.col.pos)) {
    custom.col.pos <- rep(2, length(custom.columns))
  }
  if (!is.numeric(custom.col.pos)) {
    stop("Custom column positions must be provided as a numeric vector.")
  }
  if (length(custom.col.pos) != length(custom.columns)) {
    stop(paste("Length of 'custom.col.pos' does not match length of", 
        "'custom.columns'."))
  }
  if (any(custom.col.pos > ncol(output.matrix) + 1)) {
    stop(paste("The table has only", ncol(output.matrix), "columns. The",
        "'custom.col.pos' argument does not match these dimensions."))
  }
  if (0 %in% custom.col.pos) {
    stop(paste("0 is not a valid argument in 'custom.col.pos'.", 
        "The column indices start with 1."))
  }
  for (i in 1:length(custom.columns)) {
    l <- length(custom.columns[[i]])
    if (l != numcoef && !is.null(groups) && l == (numcoef + length(groups))) {
      numcoef <- numcoef + length(groups)
    }
  }
  
  # prepare vector with column indices for custom columns
  custom.indices <- logical()
  for (i in 1:ncol(output.matrix)) {
    if (i %in% custom.col.pos) {
      custom.indices <- c(custom.indices, rep(TRUE, 
          length(which(custom.col.pos == i))), FALSE)
    } else {
      custom.indices <- c(custom.indices, FALSE)
    }
  }
  for (i in 1:length(custom.col.pos)) {
    if ((ncol(output.matrix) + 1) == custom.col.pos[i]) {
      custom.indices <- c(custom.indices, TRUE)
    }
  }
  
  if (modelnames == TRUE) {
    offset <- 1
  } else {
    offset <- 0
  }
  
  # combine output matrix with custom columns
  output.count <- 0
  custom.count <- 0
  temp <- matrix(character(), nrow = nrow(output.matrix), 
      ncol = 0)
  for (i in 1:length(custom.indices)) {
    if (custom.indices[i] == FALSE) {
      output.count <- output.count + 1
      temp <- cbind(temp, cbind(output.matrix[, output.count]))
    } else {
      custom.count <- custom.count + 1
      newcol <- matrix("", nrow = nrow(temp), ncol = 1)
      newcol[1, 1] <- names(custom.columns)[custom.count]
      for (j in 1:numcoef) {
        if (single.row == TRUE) {
          newcol[j + offset, 1] <- as.character(
              custom.columns[[custom.count]][j])
        } else {
          newcol[(2 * j) - (1 - offset), 1] <- as.character(
              custom.columns[[custom.count]][j])
        }
      }
      temp <- cbind(temp, newcol)
    }
  }
  return(temp)
}

# determine column names or column types if custom columns are present
customcolumnnames <- function(modelnames, custom.columns, custom.col.pos, 
    types = FALSE) {
  
  # adjust arguments
  modelnames <- c("", modelnames)
  if (is.null(custom.columns)) {
    if (types == FALSE) {
      return(modelnames)
    } else {
      return(c("coefnames", rep("coef", length(modelnames) - 1)))
    }
  }
  if (!class(custom.columns) == "list") {
    custom.columns <- list(custom.columns)
  }
  if (is.null(custom.col.pos) && !is.null(custom.columns)) {
    custom.col.pos <- rep(2, length(custom.columns))
  }
  
  # create indices for name injection
  custom.types <- character()
  for (i in 1:length(modelnames)) {
    if (i %in% custom.col.pos) {
      if (i == 1 && types == TRUE) {
        value <- "coefnames"
      } else {
        value <- "coef"
      }
      custom.types <- c(custom.types, rep("customcol", 
          length(which(custom.col.pos == i))), value)
    } else {
      if (i == 1) {
        custom.types <- c(custom.types, "coefnames")
      } else {
        custom.types <- c(custom.types, "coef")
      }
    }
  }
  if ((length(modelnames) + 1) %in% custom.col.pos) {
    custom.types <- c(custom.types, "customcol")
  }
  
  # do the adjustment
  output.count <- 0
  custom.count <- 0
  temp <- character()
  for (i in 1:length(custom.types)) {
    if (custom.types[i] %in% c("coef", "coefnames")) {
      output.count <- output.count + 1
      temp <- c(temp, modelnames[output.count])
    } else {
      custom.count <- custom.count + 1
      temp <- c(temp, names(custom.columns)[custom.count])
    }
  }
  
  if (types == TRUE) {
    return(custom.types)
  } else {
    return(temp)
  }
}

# print method for texreg table strings
print.texregTable <- function(x, ...) {
  cat(x, ...)
}
