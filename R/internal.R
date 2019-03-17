# The texreg package was written by Philip Leifeld.
# Please use the issue tracker at http://github.com/leifeld/texreg
# for bug reports, help or feature requests.


#' Internal functions for the \pkg{texreg} package
#'
#' @name internal
NULL


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


# function which reformats a coefficient with a certain number of decimal places
coeftostring <- function(x, lead.zero = FALSE, digits = 2) {
  if (is.na(digits)) {
    return("")
  }
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


names2latex <- function(x) {
  if (is.null(x)) {
    return(NULL)
  } else if (!grepl("\\$", x)) {
    x <- gsub("_", "\\\\_", x)
    x <- gsub("<", "\\$<\\$", x)
    x <- gsub(">", "\\$>\\$", x)
    x <- gsub("%", "\\\\%", x)
    x <- gsub("\\^2", "\\$^2\\$", x)
    x <- gsub("\\^3", "\\$^3\\$", x)
    x <- gsub("\\^4", "\\$^4\\$", x)
    x <- gsub("\\^5", "\\$^5\\$", x)
  }
  return(x)
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
            rn[i] <- gsub("\\^2", "\\$^2\\$", rn[i])
            rn[i] <- gsub("\\^3", "\\$^3\\$", rn[i])
            rn[i] <- gsub("\\^4", "\\$^4\\$", rn[i])
            rn[i] <- gsub("\\^5", "\\$^5\\$", rn[i])
        }
    }
    rownames(m) <- rn
    return(m)
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
  new.mat <- mat[new.order, , drop = FALSE]
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

# extract coefficients using the broom package
broom_coefficients <- function(x) {
    out <- broom::tidy(x)
    out <- out[, c('term', 'estimate', 'std.error', 'p.value')]
    return(out)
}

# extract gof using the broom package
broom_gof <- function(x) {
    # extract
    out <- broom::glance(x)[1, ]
    gof.decimal <- sapply(out, function(k) class(k)[1]) # type inference
    gof.decimal <- ifelse(gof.decimal %in% c('integer', 'logical'), FALSE, TRUE)
    out <- data.frame('gof.names' = colnames(out),
                      'gof' = as.numeric(out),
                      'gof.decimal' = gof.decimal,
                      stringsAsFactors = FALSE)
    # rename
    gof_dict <- c(
        'adj.r.squared' = 'Adj.\ R$^2$',
        'deviance' = 'Deviance',
        'df' = 'DF',
        'df.residual' = 'DF Resid.',
        'finTol' = 'Tolerance',
        'isConv' = 'Convergence',
        'logLik' = 'Log Likelihood',
        'null.deviance' = 'Deviance (Null)',
        'p.value' = 'P Value',
        'r.squared' = 'R$^2$',
        'sigma' = 'Sigma',
        'statistic' = 'Statistic'
    )
    gof_dict <- gof_dict[names(gof_dict) %in% out$gof.names]
    idx <- match(names(gof_dict), out$gof.names)
    out$gof.names[idx] <- gof_dict
    if (any(is.na(out$gof))) {
        warning(paste("texreg used the broom package to extract the following GOF",
                      "measures, but could not cast them to numeric type:",
                      out$gof.names[is.na(out$gof)]))
    }
    out <- stats::na.omit(out)
    return(out)
}

# create the star note (legend) printed at the bottom of tables and the stars
# printed next to standard errors
get_stars <- function(pval = NULL, # test statistics; 
                      # leave NULL if you only want the legend
                      stars = c(0.01, 0.05, 0.1),  # numeric vector of cut-offs
                      star.symbol = '*',  # character to repeat for first 3
                      # levels of significance
                      symbol = '.',  # character for 4th level of significance
                      star.prefix = '',
                      star.suffix = '',
                      ci = FALSE,  
                      ci.test = NULL, 
                      css.sup = NULL,  
                      output = 'ascii') {

    # sanity checks and prep
    if (!output %in% c('ascii', 'latex', 'html')) {
        stop("'output' argument must be 'ascii', 'latex', or 'html'.")
    }
    if (!is.numeric(ci.test) && !is.null(ci.test)) {
        stop("The argument 'ci.test' must be NULL or numeric.")
    }
    if (!is.logical(ci)) {
        stop("The argument 'ci' must be logical.")
    }
    if (!is.null(stars) && !is.numeric(stars)) {
        stop("The argument 'stars' must be NULL or numeric.")
    }
    if (any(stars > 1) | any(stars < 0)) {
        stop("All values in the 'stars' argument must be between 0 and 1.")
    }
    if (length(stars) > 4) {
        stop("Length of the 'stars' argument must be smaller than 5.") 
    }
    if (!is.null(stars) && any(is.na(stars))) {
        stop("NA value are not allowed in the 'stars' argument.")
    } 
    if (length(unique(stars)) != length(stars)) {
        stop("Duplicate elements are not allowed in the 'stars' argument.")
    }
    if (!is.null(symbol) && !is.character(symbol)) {
        stop("The argument 'symbol' must be NULL or character.")
    }
    if (is.null(css.sup) & (output == 'html')) {
        stop("To write a star note in html, you must supply 'css.sup'.")
    }
    if (length(stars) == 0) {
        stars <- NULL
    }
    p_note_flag <- any(!ci) # at least one model doesn't print confidence interval
    ci_note_flag <- any(ci) # at least one model prints a confidence interval
    
    symbols <- c(paste(rep(star.symbol, 3), collapse = ''),
                 paste(rep(star.symbol, 2), collapse = ''),
                 star.symbol,
                 symbol)
    if (length(stars) == 1) {
        symbols <- symbols[3]
    } else if (length(stars) == 2) {
        symbols <- symbols[2:3]
    } else if (length(stars) == 3) {
        symbols <- symbols[1:3]
    } 
    
    # p_note
    if (p_note_flag && !is.null(stars)) {  # stars supplied = build note
        st <- sort(stars)
        if (output == "ascii") {
            p_note <- paste0(star.prefix,
                             symbols, 
                             star.suffix,
                             " p < ", st)
        } else if (output == "latex") {
            p_note <- paste0("$^{", 
                             star.prefix,
                             symbols, 
                             star.suffix,
                             "}p<", 
                             st, 
                             "$")
        } else if (output == "html") {
            p_note <- paste0("<sup", 
                             css.sup, 
                             ">", 
                             star.prefix,
                             symbols, 
                             star.suffix,
                             "</sup>p &lt; ", 
                             st)
        }
        p_note <- paste(p_note, collapse = "; ")
    } else { # no stars supplied = empty note
        p_note <- ""
    }
    
    # ci_note
    if (ci_note_flag) {  # ci calculated for at least one model -> build ci note
        if (is.numeric(ci.test) && !is.na(ci.test)) { # sanity check  
            ci_note <- paste(ci.test, "outside the confidence interval")
        } else {
            ci_note <- ''
        }
        if (output == 'ascii') {
            ci_symbol <- '*'
        } else if (output == 'latex') {
            ci_symbol <- '$^*$'
        } else if (output == 'html') {
            ci_symbol <- paste0('<sup', css.sup, '>*</sup>') 
        }
    } else {  # ci not calculated for any model -> empty ci note
        ci_note <- ''
    }
    
    # combine p and ci notes
    if (p_note_flag && ci_note_flag && (p_note != '') && (ci_note != '')) {
        s_note <- paste0(p_note, ' (or ', ci_note, ').')
    } else if (p_note_flag && (p_note != '')) {
        s_note <- p_note
    } else if (ci_note_flag && (ci_note != '')) {
        s_note <- paste0(ci_symbol, ' ', ci_note, '.')
    } else {
        s_note <- ''
    }
    
    # stars for individual coefficients 
    if (is.null(pval) | is.null(stars)) {
        p = ''
    } else {
        if (is.na(pval)) {
            pval <- 1.0
        }
        idx <- pval < st 
        if (any(idx)) {  # choose lowest threshold (relies on previous sorting)
            p <- paste0(star.prefix, symbols[idx][1], star.suffix)
        } else {  # not significant
            p <- ""
        }
    }
    
    out <- list('note' = s_note, 'coefficients' = p)
    return(out)
}
