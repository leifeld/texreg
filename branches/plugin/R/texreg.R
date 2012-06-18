#TODO: fix all bugs; implement lm models; allow resorting the coefficient rows; implement custom coefficient labels as an argument.


# function which reformats a coefficient with two decimal places
coef.to.string <- function(x, lead.zero=FALSE) {
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


# extension for lme objects
extract.lme <- function(model) {
  
  if (!class(model) == "lme") {
    stop("Internal error: Incorrect model type! Should be an lme object!")
  }
  
  tab <- summary(model)$tTable[,-3:-4] #extract coefficient table
  
  lik <- summary(model)$logLik #extract log likelihood
  aic <- summary(model)$AIC #extract AIC
  bic <- summary(model)$BIC #extract BIC
  gof <- matrix(c(aic, bic, lik), ncol=1)
  row.names(gof) <- c("AIC", "BIC", "Log Likelihood")
  
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


texreg <- function(l, single.row=FALSE, no.margin=TRUE, leading.zero=TRUE, 
    table=TRUE, strong.signif=TRUE, symbol="\\cdot", use.packages=TRUE, 
    caption="Statistical models", label="table:coefficients", dcolumn=TRUE, 
    booktabs=TRUE, scriptsize=TRUE) {
  
  string <- ""
  
  # if a single model is handed over, put model inside a list
  if (class(l) == "ergm" | class(l) == "lme") {
    l <- list(l)
  } else if (class(l) != "list") {
    stop("Unknown object was handed over.")
  }
  
  # extract data from the models
  models <- list()
  for (i in 1:length(l)) {
    if (class(l[[i]]) == "ergm") {
      model <- extract.ergm(l[[i]])
      models <- append(models, model)
    } else if (class(l[[i]]) == "lme") {
      model <- extract.lme(l[[i]])
      models <- append(models, model)
    } else {
      warning(paste("Skipping unknown model of type ", class(l[[i]]), ".", 
          sep=""))
    }
  }
  
  # extract names of the goodness-of-fit statistics
  gof.names <- character()
  for (i in 1:length(models)) {
    for (j in 1:length(models[[i]][[2]])) {
      if (!row.names(models[[i]][[2]][j]) %in% gof.names) {
        gof.names <- append(gof.names, row.names(models[[i]][[2]][j]))
      }
    }
  }
  gof.names <- sort(gof.names)
  
  # aggregate goodness-of-fit statistics in a matrix and create list of coefs
  coefs <- list()
  gofs <- matrix(nrow=length(gof.names), ncol=length(models))
  row.names(gofs) <- gof.names
  for (i in 1:length(models)) {
    coefs <- append(coefs, models[[i]][1])
    for (j in 1:length(models[[i]][[2]])) {
      rn <- row.names(models[[i]][[2]][j])
      val <- models[[i]][[2]][j]
      col <- i
      row <- which(row.names(gofs) == rn)
      gofs[row,col] <- val
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
    if (booktabs == TRUE) {
      string <- paste(string, "\\usepackage{booktabs}\n", sep="")
    }
    if (dcolumn == TRUE) {
      string <- paste(string, "\\usepackage{dcolumn}\n\n", sep="")
    }
  }
  if (table == TRUE) {
    string <- paste(string, "\\begin{table}\n", sep="")
    string <- paste(string, "\\begin{center}\n", sep="")
    if (scriptsize == TRUE) {
      string <- paste(string, "\\scriptsize\n", sep="")
    }
  }
  string <- paste(string, "\\begin{tabular}{l ", sep="")
  
  for(i in 1:length(l)) {
    if (dcolumn == TRUE) {
      dec.left <- max(c(nchar(aics[i])-3, nchar(bics[i])-3, nchar(liks[i])-3), 3)
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
  
  if (booktabs == TRUE) {
    string <- paste(string, "}\n", "\\toprule\n", sep="")
  } else {
    string <- paste(string, "}\n", "\\hline\n", sep="")
  }

  for (k in 1:lab.length) {
    string <- paste(string, " ", sep="")
  }
  if (dcolumn == TRUE) {
    for (i in 1:length(l)) {
      string <- paste(string, " & \\multicolumn{1}{c}{Model ", i, "}", sep="")
    }
  } else {
    for (i in 1:length(l)) {
      string <- paste(string, " & Model ", i, sep="")
    }
  }
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
    for (j in 2:length(gofs[1,])) {
      string <- coef.to.string(gofs[i,j-1], leading.zero)
      gof.matrix[i,j] <- paste(dollar, string, dollar, sep="")
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
    string <- paste(string, "\\end{table}\n", sep="")
  }
  
  cat(string)
  return(string)
}
