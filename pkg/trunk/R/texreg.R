# The texreg package was written by Philip Leifeld.
# Please use the forum at http://r-forge.r-project.org/projects/texreg/ 
# for bug reports, help or feature requests.

#function which conflates a matrix with duplicate row names
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


# texreg function
texreg <- function(l, single.row=FALSE, no.margin=TRUE, leading.zero=TRUE, 
    table=TRUE, sideways=FALSE, float.pos="", strong.signif=FALSE, 
    symbol="\\cdot", use.packages=TRUE, caption="Statistical models", 
    label="table:coefficients", dcolumn=TRUE, booktabs=TRUE, scriptsize=FALSE, 
    custom.names=NA, model.names=NA, digits=2, ...) {
  
  string <- ""
  
  # if a single model is handed over, put model inside a list
  if (class(l) != "list") {
    l <- list(l)
  }

  # extract data from the models
  models <- NULL
  for (i in 1:length(l)) {
    model <- extract(l[[i]], ...)
    if (class(model[[1]]) == "list") {  #nested list of models (e.g. systemfit)
      models <- append(models, model)
    } else {  #normal case; one model
      models <- append(models, list(model))
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
  
  m <- rearrangeMatrix(m)
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
      gof.list.string[j] <- coeftostring(gof.list[j], leading.zero, 
          digits=digits)
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
    } else if (length(model.names) != length(models)) {
      stop(paste("There are", length(models), "models, but you provided", 
          length(model.names), "names for them."))
    } else {
      if (dcolumn == TRUE) {
        for (i in 1:length(models)) {
          string <- paste(string, " & \\multicolumn{1}{c}{", model.names[i], 
              "}", sep="")
        }
      } else {
        for (i in 1:length(models)) {
          string <- paste(string, " & ", model.names[i], sep="")
        }
      }
    }
  } else if (!is.na(model.names) & class(model.names) != "character") {
    stop("Model names must be specified as a vector of strings.")
  } else if (class(model.names) == "character" & 
      length(model.names) != length(models)) {
    stop(paste("A single model name was specified. But there are in fact", 
        length(models), "models."))
  } else if (class(model.names) == "character") {
    if (dcolumn == TRUE) {
      string <- paste(string, " & \\multicolumn{1}{c}{", model.names, "}", 
          sep="")
    } else {
      string <- paste(string, " & ", model.names, sep="")
    }
  } else {
    if (dcolumn == TRUE) {
      for (i in 1:length(models)) {
        string <- paste(string, " & \\multicolumn{1}{c}{Model ", i, "}", sep="")
      }
    } else {
      for (i in 1:length(models)) {
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
          output.matrix[i,k] <- "\\multicolumn{1}{c}{$-$Inf}"
        } else {
          std <- paste(" \\; (", coeftostring(m[i,j+1], leading.zero, 
              digits=digits), ")", sep="")
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
          output.matrix[(i*2)-1,k] <- "\\multicolumn{1}{c}{$-$Inf}" #upper row
          output.matrix[(i*2),k] <- "" #lower std row
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
          output.matrix[(i*2)-1,k] <- paste(dollar, coeftostring(m[i,j], 
              leading.zero, digits=digits), p, dollar, sep="")
          output.matrix[(i*2),k] <- paste(dollar, "(", coeftostring(m[i,j+1], 
              leading.zero, digits=digits), ")", dollar, sep="")
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
      strg <- coeftostring(gofs[i,j], leading.zero, digits=digits)
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
          rn == "obs" | rn == "Obs." | rn == "obs." | rn == "Events" | 
          rn == "# Events" | rn == "# events" | rn == "events" | 
          rn == "Missings" | rn == "Missing" | rn == "# Missing") {
        strg <- strsplit(strg, "\\.")[[1]][1]
        if (is.na(strg) == TRUE) {
          strg <- ""
        }
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
    string <- paste(string, "\\multicolumn{", length(models)+1, 
        "}{l}{\\textsuperscript{***}$p<0.001$, ", 
        "\\textsuperscript{**}$p<0.01$, \\textsuperscript{*}$p<0.05$, ", 
        "\\textsuperscript{$", symbol, "$}$p<0.1$}\n", sep="")
  } else {
    string <- paste(string, "\\multicolumn{", length(models)+1, 
        "}{l}{\\textsuperscript{***}$p<0.01$, ", 
        "\\textsuperscript{**}$p<0.05$, \\textsuperscript{*}$p<0.1$}\n", 
        sep="")
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
