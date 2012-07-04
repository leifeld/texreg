# The texreg package was written by Philip Leifeld.
# Please use the forum at http://r-forge.r-project.org/projects/texreg/ 
# for bug reports, help or feature requests.

# function which reformats a coefficient with two decimal places
coef.to.string <- function(x, lead.zero=FALSE) {
  if (is.na(x)) {
    return("")
  }
  y <- as.character(round(x, 2))
  if (grepl("\\.", y) == FALSE) {
    y <- paste(y, ".00", sep="")
  } else if (grepl("\\.[0-9]$", y) == TRUE) {
    y <- paste(y, "0", sep="")
  }
  if (lead.zero==FALSE && (grepl("^0",y) == TRUE || grepl("^-0",y) == TRUE)) {
    y <- gsub("0\\.", "\\.", y)
  }
  return(y)
}


# extension for lm objects
extract.lm <- function(model) {
  
  if (!class(model) == "lm") {
    stop("Internal error: Incorrect model type! Should be an lm object!")
  }
  
  tab <- summary(model)$coef[,-3] #extract coefficient table
  
  rs <- summary(model)$r.squared #extract R-squared
  adj <- summary(model)$adj.r.squared #extract adjusted R-squared
  n <- nobs(model) #extract number of observations
  
  gof <- matrix(c(rs, adj, n), ncol=1)
  row.names(gof) <- c("R$^2$", "Adj. R$^2$", "Num. obs.")
  
  table.content <- list(tab, gof)
  return(table.content)
}


# extension for gls objects
extract.gls <- function(model) {
  
  if (!class(model) == "gls") {
    stop("Internal error: Incorrect model type! Should be a gls object!")
  }
  
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


# extension for glm objects
extract.glm <- function(model) {
  
  if (!class(model)[1] == "glm") {
    stop("Internal error: Incorrect model type! Should be a glm object!")
  }
  
  tab <- summary(model)$coef[,-3] #extract coefficient table
  
  aic <- summary(model)$aic #extract AIC
  n <- nobs(model) #extract number of observations
  
  gof <- matrix(c(aic, n), ncol=1)
  row.names(gof) <- c("AIC", "Num. obs.")
  
  table.content <- list(tab, gof)
  return(table.content)
}


# extension for lme objects
extract.lme <- function(model) {
  
  if (!class(model) == "lme") {
    stop("Internal error: Incorrect model type! Should be an lme object!")
  }
  
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


# extension for ergm objects
extract.ergm <- function(model) {
  
  if (!class(model) == "ergm") {
    stop("Internal error: Incorrect model type! Should be an ergm object!")
  }
  
  tab <- summary(model)$coefs[,-3] #extract coefficient table
  
  lik <- model$mle.lik[1] #extract log likelihood
  aic <- summary(model)$aic #extract AIC
  bic <- summary(model)$bic #extract BIC
  gof <- matrix(c(aic, bic, lik), ncol=1)
  row.names(gof) <- c("AIC", "BIC", "Log Likelihood")
  
  table.content <- list(tab, gof)
  return(table.content)
}

# extension for lrm objects (Design package)
extract.lrm <- function(model) {

  if (!class(model)[1] == "lrm") {
    stop("Internal error: Incorrect model type! Should be an lrm object!")
  }

  tab <- model$coef #extract coefficient table
  tab <- cbind(COEFEST = model$coef, SE = sqrt(diag(model$var)), 
      PVALUES = pnorm(abs(model$coef/sqrt(diag(model$var))), 
      lower.tail = FALSE)*2)
  
  pseudors <- model$stats[10] #extract pseudo R-squared
  LR <- model$stats[3] #extract LR
  n <- model$stats[1] #extract number of observations
  gof <- matrix(c(pseudors, LR, n), ncol=1) #put gof measures in a 1-column matrix
  row.names(gof) <- c("Pseudo R$^2$", "L.R.", "Num. obs.") #set row names
  table.content <- list(tab, gof) #put coefficients and gofs in a list
  
  return(table.content) #return the list object
}

texreg <- function(l, single.row=FALSE, no.margin=TRUE, leading.zero=TRUE, 
    table=TRUE, sideways=FALSE, float.pos="", strong.signif=FALSE, 
    symbol="\\cdot", use.packages=TRUE, caption="Statistical models", 
    label="table:coefficients", dcolumn=TRUE, booktabs=TRUE, scriptsize=FALSE, 
    custom.names=NA, model.names=NA) {
  
  string <- ""
  
  # if a single model is handed over, put model inside a list
  # IMPLEMENT NEW EXTENSIONS HERE
  if (class(l)[1] == "ergm" | class(l)[1] == "lme" | class(l)[1] == "lm" | 
      class(l)[1] == "gls" | class(l)[1] == "glm" | class(l)[1] == "lrm") {
    l <- list(l)
  } else if (class(l) != "list") {
    stop("Unknown object was handed over.")
  }
  
  # extract data from the models
  # IMPLEMENT NEW EXTENSIONS HERE
  models <- NULL
  for (i in 1:length(l)) {
    if (class(l[[i]])[1] == "ergm") {
      model <- extract.ergm(l[[i]])
      models <- append(models, list(model))
    } else if (class(l[[i]])[1] == "lme") {
      model <- extract.lme(l[[i]])
      models <- append(models, list(model))
    } else if (class(l[[i]])[1] == "lm") {
      model <- extract.lm(l[[i]])
      models <- append(models, list(model))
    } else if (class(l[[i]])[1] == "gls") {
      model <- extract.gls(l[[i]])
      models <- append(models, list(model))
    } else if (class(l[[i]])[1] == "glm") {
      model <- extract.glm(l[[i]])
      models <- append(models, list(model))
    } else if (class(l[[i]])[1] == "lrm") {
      model <- extract.lrm(l[[i]])
      models <- append(models, list(model))
    } else {
      warning(paste("Skipping unknown model of type ", class(l[[i]]), ".", 
          sep=""))
    }
  }
  
  # extract names of the goodness-of-fit statistics
  gof.names <- character()
  for (i in 1:length(models)) {
    for (j in 1:length(models[[i]][[2]])) {
      if (!row.names(models[[i]][[2]])[j] %in% gof.names) {
        gof.names <- append(gof.names, row.names(models[[i]][[2]])[j])
      }
    }
  }
  
  # aggregate goodness-of-fit statistics in a matrix and create list of coefs
  coefs <- list()
  gofs <- matrix(nrow=length(gof.names), ncol=length(models))
  row.names(gofs) <- gof.names
  for (i in 1:length(models)) {
    coefs <- append(coefs, models[[i]][1])
    for (j in 1:length(models[[i]][[2]])) {
      rn <- row.names(models[[i]][[2]])[j]
      val <- models[[i]][[2]][[j]]
      col <- i
      row <- which(row.names(gofs) == rn)
      gofs[row,col] <- val
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
  
  # use custom coefficient names if provided
  if (length(custom.names) > 1) {
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
  
  # check if the custom name procedure caused duplicates and merge them
  for (i in 1:length(rownames(m))) {  #go through rows
    for (j in 1:length(rownames(m))) {  #go through rows (to find duplicates)
      if (i != j & rownames(m)[i] == rownames(m)[j]) {  #found a duplicate name
        identical <- logical(length(m[i,]))
        for (k in 1:length(m[i,])) {  #go through columns
          if ( (is.na(m[i,k]) & !is.na(m[j,k])) | 
              (!is.na(m[i,k]) & is.na(m[j,k])) | 
              (is.na(m[i,k]) & is.na(m[j,k])) ) {
            identical[k] <- TRUE  #set TRUE if they are complementary
          }
        }
        if (length(identical[identical==FALSE]) == 0) {  #if complementary...
          for (k in 1:ncol(m)) {  #go through the columns again
            if (is.na(m[i,k])) {
              m[i,k] <- m[j,k]  #merge them
            } else if (is.na(m[j,k])) {
              m[j,k] <- m[i,k]  #merge them
            }
          }
        }
      }
    }
  }
  m <- m[duplicated(m) == FALSE,]  #remove duplicate rows
  m <- as.data.frame(m)
  
  # what is the optimal length of the labels?
  lab.list <- c(rownames(m), gof.names)
  lab.length <- 0
  for (i in 1:length(lab.list)) {
    if (nchar(lab.list[i]) > lab.length) {
      lab.length <- nchar(lab.list[i])
    }
  }
  
  # write table header
  string <- paste(string, "\n", sep="")
  if (use.packages == TRUE) {
    if (sideways == TRUE & table == TRUE) {
      string <- paste(string, "\\usepackage{rotating}\n", sep="")
    }
    if (booktabs == TRUE) {
      string <- paste(string, "\\usepackage{booktabs}\n", sep="")
    }
    if (dcolumn == TRUE) {
      string <- paste(string, "\\usepackage{dcolumn}\n\n", sep="")
    }
  }
  if (table == TRUE) {
    if (sideways == TRUE) {
      t <- "sideways"
    } else {
      t <- ""
    }
    if ( float.pos == "") {
      string <- paste(string, "\\begin{", t, "table}\n", sep="")
    } else {
      string <- paste(string, "\\begin{", t, "table}[", float.pos, "]\n", 
          sep="")
    }
    string <- paste(string, "\\begin{center}\n", sep="")
    if (scriptsize == TRUE) {
      string <- paste(string, "\\scriptsize\n", sep="")
    }
  }
  string <- paste(string, "\\begin{tabular}{l ", sep="")
  
  #define columns of the table
  for(i in 1:length(models)) {
    gof.list <- as.vector(gofs[,i])
    gof.list.string <- NULL
    for (j in 1:length(gof.list)) {
      gof.list.string[j] <- coef.to.string(gof.list[j], leading.zero)
    }
    if (dcolumn == TRUE) {
      dec.left <- max(c(nchar(gof.list.string)-3), 3)
      if (single.row == TRUE) {
        dec.right <- 3
        separator <- ")"
        dec.left <- 11
      } else {
        dec.right <- 5
        separator <- "."
      }
      if (no.margin == FALSE) {
        margin.arg <- ""
      } else {
        margin.arg <- "@{}"
      }
      string <- paste(string, "D{", separator, "}{", separator, "}{", dec.left, 
          separator, dec.right, "} ", margin.arg, sep="")
    } else {
      string <- paste(string, "c ", sep="")
    }
  }
  
  # horizontal rule above the table
  if (booktabs == TRUE) {
    string <- paste(string, "}\n", "\\toprule\n", sep="")
  } else {
    string <- paste(string, "}\n", "\\hline\n", sep="")
  }
  
  # specify model names
  for (k in 1:lab.length) {
    string <- paste(string, " ", sep="")
  }
  if (length(model.names) > 1) {
    if (class(model.names) != "character") {
      stop("Model names must be specified as a vector of strings.")
    } else if (length(model.names) != length(l)) {
      stop(paste("There are", length(l), "models, but you provided", 
          length(model.names), "names for them."))
    } else {
      if (dcolumn == TRUE) {
        for (i in 1:length(l)) {
          string <- paste(string, " & \\multicolumn{1}{c}{", model.names[i], 
              "}", sep="")
        }
      } else {
        for (i in 1:length(l)) {
          string <- paste(string, " & ", model.names[i], sep="")
        }
      }
    }
  } else if (!is.na(model.names) & class(model.names) != "character") {
    stop("Model names must be specified as a vector of strings.")
  } else if (class(model.names) == "character" & 
      length(model.names) != length(l)) {
    stop(paste("A single model name was specified. But there are in fact", 
        length(l), "models."))
  } else if (class(model.names) == "character") {
    if (dcolumn == TRUE) {
      string <- paste(string, " & \\multicolumn{1}{c}{", model.names, "}", 
          sep="")
    } else {
      string <- paste(string, " & ", model.names, sep="")
    }
  } else {
    if (dcolumn == TRUE) {
      for (i in 1:length(l)) {
        string <- paste(string, " & \\multicolumn{1}{c}{Model ", i, "}", sep="")
      }
    } else {
      for (i in 1:length(l)) {
        string <- paste(string, " & Model ", i, sep="")
      }
    }
  }
  
  # horizontal rule between coefficients and goodness-of-fit block
  if (booktabs == TRUE) {
    string <- paste(string, " \\\\\n", "\\midrule\n", sep="")
  } else {
    string <- paste(string, " \\\\\n", "\\hline\n", sep="")
  }
  
  # write coefficient rows
  if (single.row==TRUE) {
    output.matrix <- matrix(ncol=(length(m)/3)+1, nrow=length(m[,1]))
    
    # row labels
    for (i in 1:length(rownames(m))) {
      output.matrix[i,1] <- rownames(m)[i]
    }
    
    # coefficients and standard deviations
    for (i in 1:length(m[,1])) { #go through rows
      j <- 1 #column in the original, merged coef table
      k <- 2 #second column of output.matrix, i.e., coefficients
      while (j <= length(m)) {
        if (is.na(m[i,j])) {
          output.matrix[i,k] <- ""
        } else if (m[i,j] == -Inf) {
          output.matrix[i,k] <- "-Inf (NA)"
        } else {
          std <- paste(" \\; (", coef.to.string(m[i,j+1], leading.zero), ")", 
              sep="")
          if (strong.signif == TRUE) {
            if (m[i,j+2] <= 0.001) {
              p <- "^{***}"
            } else if (m[i,j+2] <= 0.01) {
              p <- "^{**}"
            } else if (m[i,j+2] <= 0.05) {
              p <- "^{*}"
            } else if (m[i,j+2] <= 0.1) {
              p <- paste("^{", symbol, "}", sep="")
            } else {
              p <- ""
            }
          } else {
            if (m[i,j+2] <= 0.01) {
              p <- "^{***}"
            } else if (m[i,j+2] <= 0.05) {
              p <- "^{**}"
            } else if (m[i,j+2] <= 0.1) {
              p <- "^{*}"
            } else {
              p <- ""
            }
          }
          if (dcolumn == TRUE) {
            dollar <- ""
          } else {
            dollar <- "$"
          }
          entry <- paste(dollar, coef.to.string(m[i,j], leading.zero), std, p, 
              dollar, sep="")
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
        if (is.na(m[i,j])) {
          output.matrix[(i*2)-1,k] <- "" #upper coefficient row
          output.matrix[(i*2),k] <- "" #lower std row
        } else if (m[i,j] == -Inf) {
          output.matrix[(i*2)-1,k] <- "-Inf" #upper coefficient row
          output.matrix[(i*2),k] <- "(NA)" #lower std row
        } else {
          if (strong.signif == TRUE) {
            if (m[i,j+2] <= 0.001) {
              p <- "^{***}"
            } else if (m[i,j+2] <= 0.01) {
              p <- "^{**}"
            } else if (m[i,j+2] <= 0.05) {
              p <- "^{*}"
            } else if (m[i,j+2] <= 0.1) {
              p <- paste("^{", symbol, "}", sep="")
            } else {
              p <- ""
            }
          } else {
            if (m[i,j+2] <= 0.01) {
              p <- "^{***}"
            } else if (m[i,j+2] <= 0.05) {
              p <- "^{**}"
            } else if (m[i,j+2] <= 0.1) {
              p <- "^{*}"
            } else {
              p <- ""
            }
          }
          if (dcolumn == TRUE) {
            dollar <- ""
          } else {
            dollar <- "$"
          }
          output.matrix[(i*2)-1,k] <- paste(dollar, coef.to.string(m[i,j], 
              leading.zero), p, dollar, sep="")
          output.matrix[(i*2),k] <- paste(dollar, "(", coef.to.string(m[i,j+1], 
              leading.zero), ")", dollar, sep="")
        }
        k <- k+1
        j <- j+3
      }
    }
  }
  
  # goodness-of-fit statistics
  if (dcolumn == TRUE) {
    dollar <- ""
  } else {
    dollar <- "$"
  }
  gof.matrix <- matrix(nrow=nrow(gofs), ncol=ncol(gofs)+1) #incl. labels
  for (i in 1:length(gofs[,1])) {
    gof.matrix[i,1] <- rownames(gofs)[i]
    for (j in 1:length(gofs[1,])) {
      strg <- coef.to.string(gofs[i,j], leading.zero)
      rn <- rownames(gofs)[i]
      if (rn == "Num. obs." | rn == "n" | rn == "N" | rn == "N obs" | 
          rn == "N obs." | rn == "nobs" | rn == "n obs" | rn == "n obs." | 
          rn == "n.obs." | rn == "N.obs." | rn == "N. obs" | 
          rn == "Num observations" | rn == "Number of observations" | 
          rn == "Num obs" | rn == "num obs" | rn == "Num. observations" | 
          rn == "Num Observations" | rn == "Num. Observations" | 
          rn == "Num. Obs." | rn == "Num.Obs." | rn == "Number obs." | 
          rn == "Number Obs." | rn == "Number obs" | rn == "Number Obs" | 
          rn == "Number of Obs." | rn == "Number of obs." | 
          rn == "Number of obs" | rn == "Number of Obs" | rn == "Obs" | 
          rn == "obs" | rn == "Obs." | rn == "obs.") {
        strg <- substring(strg, 1, nchar(strg)-3)
      }
      gof.matrix[i,j+1] <- paste(dollar, strg, dollar, sep="")
    }
  }
  
  # combine the coefficient and gof matrices vertically
  output.matrix <- rbind(output.matrix, gof.matrix)
  
  # fill with spaces
  max.lengths <- numeric(length(output.matrix[1,]))
  for (i in 1:length(output.matrix[1,])) {
    max.length <- 0
    for (j in 1:length(output.matrix[,1])) {
      if (nchar(output.matrix[j,i]) > max.length) {
        max.length <- nchar(output.matrix[j,i])
      }
    }
    max.lengths[i] <- max.length
  }
  for (i in 1:length(output.matrix[,1])) {
    for (j in 1:length(output.matrix[1,])) {
      nzero <- max.lengths[j] - nchar(output.matrix[i,j])
      zeros <- rep(" ", nzero)
      zeros <- paste(zeros, collapse="")
      output.matrix[i,j] <- paste(output.matrix[i,j], zeros, sep="")
    }
  }
  
  
  # write coefficients to string object
  for (i in 1:(length(output.matrix[,1])-length(gof.names))) {
    for (j in 1:length(output.matrix[1,])) {
      string <- paste(string, output.matrix[i,j], sep="")
      if (j == length(output.matrix[1,])) {
        string <- paste(string, " \\\\\n", sep="")
      } else {
        string <- paste(string, " & ", sep="")
      }
    }
  }
  
  if (booktabs == TRUE) {
    string <- paste(string, "\\midrule\n", sep="")
  } else {
    string <- paste(string, "\\hline\n", sep="")
  }
  
  for (i in (length(output.matrix[,1])-(length(gof.names)-1)):
      (length(output.matrix[,1]))) {
    for (j in 1:length(output.matrix[1,])) {
      string <- paste(string, output.matrix[i,j], sep="")
      if (j == length(output.matrix[1,])) {
        string <- paste(string, " \\\\\n", sep="")
      } else {
        string <- paste(string, " & ", sep="")
      }
    }
  }
  
  # write table footer
  if (booktabs == TRUE) {
    string <- paste(string, "\\bottomrule\n", sep="")
  } else {
    string <- paste(string, "\\hline\n", sep="")
  }
  string <- paste(string, "\\vspace{-2mm}\\\\\n", sep="")
  
  if (strong.signif == TRUE) {
    string <- paste(string, "\\multicolumn{", length(l)+1, 
        "}{l}{\\textsuperscript{***}$p<0.001$, ", 
        "\\textsuperscript{**}$p<0.01$, \\textsuperscript{*}$p<0.05$, ", 
        "\\textsuperscript{$", symbol, "$}$p<0.1$}\n", sep="")
  } else {
    string <- paste(string, "\\multicolumn{", length(l)+1, 
        "}{l}{\\textsuperscript{***}$p<0.01$, ", 
        "\\textsuperscript{**}$p<0.05$, \\textsuperscript{*}$p<0.1$}\n", sep="")
  }
  
  string <- paste(string, "\\end{tabular}\n", sep="")
  
  if (table == TRUE) {
    if (scriptsize == TRUE) {
      string <- paste(string, "\\normalsize\n", sep="")
    }
    string <- paste(string, "\\end{center}\n", sep="")
    string <- paste(string, "\\caption{", caption, "}\n", sep="")
    string <- paste(string, "\\label{", label, "}\n", sep="")
    if (sideways == TRUE) {
      t <- "sideways"
    } else {
      t <- ""
    }
    string <- paste(string, "\\end{", t, "table}\n", sep="")
  }
  
  cat(string)
  return(string)
}
