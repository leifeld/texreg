# The texreg package was written by Philip Leifeld.
# Please use the forum at http://r-forge.r-project.org/projects/texreg/ 
# for bug reports, help or feature requests.


# screenreg function
screenreg <- function(l, single.row=FALSE, leading.zero=TRUE, 
    stars=c(0.001, 0.01, 0.05, 0.1), symbol=".", custom.names=NULL, 
    custom.gof.names=NULL, model.names=NULL, digits=2, outer.rule="=", 
    inner.rule="-", column.spacing=2, override.coef=0, override.se=0, 
    override.pval=0, omit.coef=NA, reorder.coef=NULL, reorder.gof=NULL, 
    file=NA, return.string=FALSE, custom.note=NULL, ...) {
  
  stars <- check.stars(stars)
  
  models <- get.data(l, ...) #extract relevant coefficients, SEs, GOFs, etc.
  models <- override(models, override.coef, override.se, override.pval)
  models <- tex.replace(models, type="screen") #convert TeX code to text code
  gof.names <- get.gof(models) #extract names of GOFs
  
  # arrange coefficients and GOFs nicely in a matrix
  gofs <- aggregate.matrix(models, gof.names, custom.gof.names, digits, 
      returnobject="gofs")
  m <- aggregate.matrix(models, gof.names, custom.gof.names, digits, 
      returnobject="m")
  decimal.matrix <- aggregate.matrix(models, gof.names, custom.gof.names, 
      digits, returnobject="decimal.matrix")
  
  m <- customnames(m, custom.names) #rename coefficients
  m <- rearrangeMatrix(m) #resort matrix and conflate duplicate entries
  m <- as.data.frame(m)
  m <- omitcoef(m, omit.coef) #remove coefficient rows matching regex
  
  modnames <- modelnames(models, model.names) #use (custom) model names
  
  # reorder GOF and coef matrix
  m <- reorder(m, reorder.coef)
  gofs <- reorder(gofs, reorder.gof)
  decimal.matrix <- reorder(decimal.matrix, reorder.gof)
  
  # create output table with significance stars etc.
  output.matrix <- outputmatrix(m, single.row, neginfstring="-Inf", 
      leading.zero, digits, se.prefix=" (", se.suffix=")", star.prefix=" ", 
      star.suffix="", star.char="*", stars, dcolumn=TRUE, symbol=symbol, 
      bold=0, bold.prefix="", bold.suffix="")
  
  # create GOF matrix (the lower part of the final output matrix)
  gof.matrix <- gofmatrix(gofs, decimal.matrix, dcolumn=TRUE, leading.zero, 
      digits)
  
  # combine the coefficient and gof matrices vertically
  output.matrix <- rbind(output.matrix, gof.matrix)
  
  # reformat output matrix and add spaces
  if (ncol(output.matrix) == 2) {
    temp <- matrix(format.column(output.matrix[,-1], single.row=single.row, 
        digits=digits))
  } else {
    temp <- apply(output.matrix[,-1], 2, format.column, single.row=single.row, 
        digits=digits)
  }
  output.matrix <- cbind(output.matrix[,1], temp)
  output.matrix <- rbind(c("", modnames), output.matrix)
  for (i in 1:ncol(output.matrix)) {
    output.matrix[,i] <- fill.spaces(output.matrix[,i])
  }
  
  string <- "\n"
  
  # horizontal rule above the table
  table.width <- sum(nchar(output.matrix[1,])) + 
      (ncol(output.matrix) - 1) * column.spacing
  if (class(outer.rule) != "character") {
    stop("outer.rule must be a character.")
  } else if (nchar(outer.rule) > 1) {
    stop("outer.rule must be a character of maximum length 1.")
  } else if (outer.rule == "") {
    o.rule <- ""
  } else {
    o.rule <- paste(rep(outer.rule, table.width), collapse="")
    string <- paste(string, o.rule, "\n", sep="")
  }
  
  # specify model names
  spacing <- paste(rep(" ", column.spacing), collapse="")
  string <- paste(string, output.matrix[1,1], sep="")
  for (i in 2:ncol(output.matrix)) {
    string <- paste(string, spacing, output.matrix[1,i], sep="")
  }
  string <- paste(string, "\n", sep="")
  
  # mid rule 1
  if (class(inner.rule) != "character") {
    stop("inner.rule must be a character.")
  } else if (nchar(inner.rule) > 1) {
    stop("inner.rule must be a character of maximum length 1.")
  } else if (inner.rule == "") {
    i.rule <- ""
  } else {
    i.rule <- paste(rep(inner.rule, table.width), collapse="")
    string <- paste(string, i.rule, "\n", sep="")
  }
  
  # write coefficients
  for (i in 2:(length(output.matrix[,1])-length(gof.names))) {
    for (j in 1:length(output.matrix[1,])) {
      string <- paste(string, output.matrix[i,j], sep="")
      if (j == length(output.matrix[1,])) {
        string <- paste(string, "\n", sep="")
      } else {
        string <- paste(string, spacing, sep="")
      }
    }
  }
  
  # mid rule 2
  if (inner.rule != "") {
    string <- paste(string, i.rule, "\n", sep="")
  }
  
  # write GOF part of the output matrix
  for (i in (length(output.matrix[,1])-(length(gof.names)-1)):
      (length(output.matrix[,1]))) {
    for (j in 1:length(output.matrix[1,])) {
      string <- paste(string, output.matrix[i,j], sep="")
      if (j == length(output.matrix[1,])) {
        string <- paste(string, "\n", sep="")
      } else {
        string <- paste(string, spacing, sep="")
      }
    }
  }
  
  # write table footer
  if (outer.rule != "") {
    string <- paste(string, o.rule, "\n", sep="")
  }
  
  # stars note
  if (!is.null(custom.note)) {
    if (custom.note == "") {
      note <- "\n"
    } else {
      note <- paste(custom.note, "\n\n", sep="")
    }
  } else if (is.null(stars)) {
    note <- "\n"
  } else {
    st <- sort(stars)
    if (length(unique(st)) != length(st)) {
      stop("Duplicate elements are not allowed in the stars argument.")
    }
    if (length(st) == 4) {
      note <- paste("*** p < ", st[1], ", ** p < ", st[2], ", * p < ", st[3], 
          ", ", symbol, " p < ", st[4], "\n\n", sep="")
    } else if (length(st) == 3) {
      note <- paste("*** p < ", st[1], ", ** p < ", st[2], ", * p < ", st[3], 
          "\n\n", sep="")
    } else if (length(st) == 2) {
      note <- paste("** p < ", st[1], ", * p < ", st[2], "\n\n", sep="")
    } else if (length(st) == 1) {
      note <- paste("* p < ", st, "\n\n", sep="")
    } else {
      note <- "\n"
    }
  }
  string <- paste(string, note, sep="")
  
  #write to file
  if (is.na(file)) {
    cat(string)
  } else if (!is.character(file)) {
    stop("The 'file' argument must be a character string.")
  } else {
    sink(file)
    cat(string)
    sink()
    cat(paste("The table was written to the file '", file, "'.\n", sep=""))
  }
  
  if (return.string == TRUE) {
    return(string)
  }
}


# texreg function
texreg <- function(l, single.row=FALSE, no.margin=TRUE, leading.zero=TRUE, 
    table=TRUE, sideways=FALSE, float.pos="", stars=c(0.001, 0.01, 0.05, 0.1), 
    symbol="\\cdot", use.packages=TRUE, caption="Statistical models", 
    label="table:coefficients", dcolumn=TRUE, booktabs=TRUE, scriptsize=FALSE, 
    custom.names=NULL, custom.gof.names=NULL, model.names=NULL, digits=2, 
    center=TRUE, override.coef=0, override.se=0, override.pval=0, omit.coef=NA, 
    reorder.coef=NULL, reorder.gof=NULL, file=NA, return.string=FALSE, 
    caption.above=FALSE, bold=0.00, custom.note=NULL, ...) {
  
  stars <- check.stars(stars)
  
  #check dcolumn vs. bold
  if (dcolumn==TRUE && bold > 0) {
    dcolumn <- FALSE
    msg <- paste("The dcolumn package and the bold argument cannot be used at", 
        "the same time. Switching off dcolumn.")
    if (stars == TRUE) {
      warning(paste(msg, "You should also consider setting stars=FALSE."))
    } else {
      warning(msg)
    }
  }
  
  models <- get.data(l, ...) #extract relevant coefficients, SEs, GOFs, etc.
  gof.names <- get.gof(models) #extract names of GOFs
  models <- override(models, override.coef, override.se, override.pval)
  
  # arrange coefficients and GOFs nicely in a matrix
  gofs <- aggregate.matrix(models, gof.names, custom.gof.names, digits, 
      returnobject="gofs")
  m <- aggregate.matrix(models, gof.names, custom.gof.names, digits, 
      returnobject="m")
  decimal.matrix <- aggregate.matrix(models, gof.names, custom.gof.names, 
      digits, returnobject="decimal.matrix")
  
  m <- customnames(m, custom.names) #rename coefficients
  m <- rearrangeMatrix(m) #resort matrix and conflate duplicate entries
  m <- as.data.frame(m)
  m <- omitcoef(m, omit.coef) #remove coefficient rows matching regex
  
  modnames <- modelnames(models, model.names) #use (custom) model names
  
  # reorder GOF and coef matrix
  m <- reorder(m, reorder.coef)
  gofs <- reorder(gofs, reorder.gof)
  decimal.matrix <- reorder(decimal.matrix, reorder.gof)
  
  # what is the optimal length of the labels?
  lab.list <- c(rownames(m), gof.names)
  lab.length <- 0
  for (i in 1:length(lab.list)) {
    if (nchar(lab.list[i]) > lab.length) {
      lab.length <- nchar(lab.list[i])
    }
  }
  
  string <- ""
  
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
    if (caption.above == TRUE) {
      string <- paste(string, "\\caption{", caption, "}\n", sep="")
    }
    if (center == TRUE) {
      string <- paste(string, "\\begin{center}\n", sep="")
    }
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

  if (dcolumn == TRUE) {
    for (i in 1:length(models)) {
      string <- paste(string, " & \\multicolumn{1}{c}{", modnames[i], 
          "}", sep="")
    }
  } else {
    for (i in 1:length(models)) {
      string <- paste(string, " & ", modnames[i], sep="")
    }
  }
  
  # horizontal rule between coefficients and goodness-of-fit block
  if (booktabs == TRUE) {
    string <- paste(string, " \\\\\n", "\\midrule\n", sep="")
  } else {
    string <- paste(string, " \\\\\n", "\\hline\n", sep="")
  }
  
  # create output table with significance stars etc.
  output.matrix <- outputmatrix(m, single.row, 
      neginfstring="\\multicolumn{1}{c}{$-$Inf}", leading.zero, digits, 
      se.prefix=" \\; (", se.suffix=")", star.prefix="^{", star.suffix="}", 
      star.char="*", stars, dcolumn=dcolumn, symbol, bold, 
      bold.prefix="\\textbf{", bold.suffix="}")
  
  # create GOF matrix (the lower part of the final output matrix)
  gof.matrix <- gofmatrix(gofs, decimal.matrix, dcolumn=TRUE, leading.zero, 
      digits)
  
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
  if (is.null(custom.note) || custom.note != "") {
    string <- paste(string, "\\vspace{-3mm}\\\\\n", sep="")
  }
  
  # stars note
  if (!is.null(custom.note)) {
    if (custom.note == "") {
      note <- ""
    } else {
      note <- paste(custom.note, "\n", sep="")
    }
  } else if (is.null(stars)) {
    note <- "\n"
  } else {
    st <- sort(stars)
    if (length(unique(st)) != length(st)) {
      stop("Duplicate elements are not allowed in the stars argument.")
    }
    if (length(st) == 4) {
      note <- paste("\\multicolumn{", length(models)+1, 
        "}{l}{\\textsuperscript{***}$p<", st[1], 
        "$, \n  \\textsuperscript{**}$p<", st[2], 
        "$, \n  \\textsuperscript{*}$p<", st[3], 
        "$, \n  \\textsuperscript{$", symbol, "$}$p<", st[4], "$}\n", sep="")
    } else if (length(st) == 3) {
      note <- paste("\\multicolumn{", length(models)+1, 
        "}{l}{\\textsuperscript{***}$p<", st[1], 
        "$, \n  \\textsuperscript{**}$p<", st[2], 
        "$, \n  \\textsuperscript{*}$p<", st[3], "$}\n", sep="")
    } else if (length(st) == 2) {
      note <- paste("\\multicolumn{", length(models)+1, 
        "}{l}{\\textsuperscript{**}$p<", st[1], 
        "$, \n  \\textsuperscript{*}$p<", st[2], "$}\n", sep="")
    } else if (length(st) == 1) {
      note <- paste("\\multicolumn{", length(models)+1, 
        "}{l}{\\textsuperscript{*}$p<", st[1], "$}\n", sep="")
    } else {
      note <- ""
    }
  }
  string <- paste(string, note, sep="")
  
  string <- paste(string, "\\end{tabular}\n", sep="")
  
  if (table == TRUE) {
    if (scriptsize == TRUE) {
      string <- paste(string, "\\normalsize\n", sep="")
    }
    if (caption.above == FALSE) {
      string <- paste(string, "\\caption{", caption, "}\n", sep="")
    }
    string <- paste(string, "\\label{", label, "}\n", sep="")
    if (center == TRUE) {
      string <- paste(string, "\\end{center}\n", sep="")
    }
    if (sideways == TRUE) {
      t <- "sideways"
    } else {
      t <- ""
    }
    string <- paste(string, "\\end{", t, "table}\n\n", sep="")
  }
  
  if (is.na(file)) {
    cat(string)
  } else if (!is.character(file)) {
    stop("The 'file' argument must be a character string.")
  } else {
    sink(file)
    cat(string)
    sink()
    cat(paste("The table was written to the file '", file, "'.\n", sep=""))
  }
  if (return.string == TRUE) {
    return(string)
  }
}


# htmlreg function
htmlreg <- function(l, single.row=FALSE, leading.zero=TRUE, 
    stars=c(0.001, 0.01, 0.05, 0.1), symbol="&middot;", caption="", 
    custom.names=NULL, custom.gof.names=NULL, model.names=NULL, digits=2, 
    doctype=TRUE, star.symbol="*", center=FALSE, override.coef=0, 
    override.se=0, override.pval=0, omit.coef=NA, reorder.coef=NULL, 
    reorder.gof=NULL, file=NA, return.string=FALSE, caption.above=FALSE, 
    bold=0.00, custom.note=NULL, ...) {
  
  stars <- check.stars(stars)
  
  models <- get.data(l, ...) #extract relevant coefficients, SEs, GOFs, etc.
  
  models <- override(models, override.coef, override.se, override.pval)
  models <- tex.replace(models, type="html") #convert TeX code to HTML code
  gof.names <- get.gof(models) #extract names of GOFs
  
  # arrange coefficients and GOFs nicely in a matrix
  gofs <- aggregate.matrix(models, gof.names, custom.gof.names, digits, 
      returnobject="gofs")
  m <- aggregate.matrix(models, gof.names, custom.gof.names, digits, 
      returnobject="m")
  decimal.matrix <- aggregate.matrix(models, gof.names, custom.gof.names, 
      digits, returnobject="decimal.matrix")
  
  m <- customnames(m, custom.names) #rename coefficients
  m <- rearrangeMatrix(m) #resort matrix and conflate duplicate entries
  m <- as.data.frame(m)
  m <- omitcoef(m, omit.coef) #remove coefficient rows matching regex
  
  modnames <- modelnames(models, model.names) #use (custom) model names
  
  # reorder GOF and coef matrix
  m <- reorder(m, reorder.coef)
  gofs <- reorder(gofs, reorder.gof)
  decimal.matrix <- reorder(decimal.matrix, reorder.gof)
  
  # create output table with significance stars etc.
  output.matrix <- outputmatrix(m, single.row, neginfstring="-Inf", 
      leading.zero, digits, se.prefix=" (", se.suffix=")", 
      star.char=star.symbol, star.prefix="<sup>", star.suffix="</sup>", 
      stars, dcolumn=TRUE, symbol, bold, bold.prefix="<b>", bold.suffix="</b>")
  
  # create GOF matrix (the lower part of the final output matrix)
  gof.matrix <- gofmatrix(gofs, decimal.matrix, leading.zero, 
      digits)
  
  # combine the coefficient and gof matrices vertically
  output.matrix <- rbind(output.matrix, gof.matrix)
  
  # write table header
  if (single.row == TRUE) {
    numcols <- 2 * length(models)
  } else {
    numcols <- length(models)
  }
  
  if (doctype == TRUE) {
    doct <- "<!DOCTYPE html>\n"
  } else {
    doct <- ""
  }
  
  if (center == FALSE) {
    tabdef <- "    <table cellspacing=\"0\">\n"
  } else {
    tabdef <- "    <table cellspacing=\"0\" align=\"center\">\n"
  }
  
  if (is.null(caption) || !is.character(caption)) {
    stop("The caption must be provided as a (possibly empty) character vector.")
  } else if (caption != "" && caption.above == FALSE) {
    cap <- paste("      <caption align=\"bottom\" style=\"margin-top:0.3em\">", 
        caption, "</caption>\n", sep="")
  } else if (caption != "" && caption.above == TRUE) {
    cap <- paste("      <caption align=\"top\" style=\"margin-bottom:0.3em\">", 
        caption, "</caption>\n", sep="")
  } else {
    cap <- ""
  }
  
  string <- paste(
      "\n", 
      doct, 
      "<html>\n", 
      "  <head>\n", 
      "    <title>", caption, "</title>\n", 
      "    <style type=\"text/css\">\n", 
      "      table {\n", 
      "        border: none;\n", 
      "      }\n", 
      "      th {\n", 
      "        text-align: left;\n", 
      "        border-top: 2px solid black;\n", 
      "        border-bottom: 1px solid black;\n", 
      "        padding-right: 12px;\n", 
      "      }\n", 
      "      .midRule {\n", 
      "        border-top: 1px solid black;\n", 
      "      }\n", 
      "      .bottomRule {\n", 
      "        border-bottom: 2px solid black;\n", 
      "      }\n", 
      "      td {\n", 
      "        padding-right: 12px;\n", 
      "        border: none;\n", 
      "      }\n", 
      "      caption span {\n", 
      "        float: left;\n", 
      "      }\n", 
      "      sup {\n", 
      "        vertical-align: 4px;\n", 
      "      }\n", 
      "    </style>\n", 
      "  </head>\n\n", 
      "  <body>\n", 
      tabdef, 
      cap, 
      "      <tr>\n", 
      "        <th class=\"modelnames\"></th>\n", 
      sep="")
  
  # specify model names (header row)
  for (i in 1:length(models)) {
    string <- paste(string, 
        "        <th><b>", 
        modnames[i], "</b></th>\n", sep="")
  }
  string <- paste(string, "      </tr>\n", sep="")
  
  # write coefficients to string object
  for (i in 1:(length(output.matrix[,1])-length(gof.names))) {
    string <- paste(string, "      <tr>\n", sep="")
    for (j in 1:length(output.matrix[1,])) {
      string <- paste(string, "        <td>", 
          output.matrix[i,j], "</td>\n", sep="")
    }
    string <- paste(string, "      </tr>\n", sep="")
  }
  
  for (i in (length(output.matrix[,1])-(length(gof.names)-1)):
      (length(output.matrix[,1]))) {
    string <- paste(string, "      <tr>\n", sep="")
    for (j in 1:length(output.matrix[1,])) {
      if (i == length(output.matrix[,1])-(length(gof.names)-1)) {
        string <- paste(string, "        <td class=\"midRule\">", 
            output.matrix[i,j], "</td>\n", sep="")
      } else if (i == length(output.matrix[,1])) {
        string <- paste(string, "        <td class=\"bottomRule\">", 
            output.matrix[i,j], "</td>\n", sep="")
      } else {
        string <- paste(string, "        <td>", 
            output.matrix[i,j], "</td>\n", sep="")
      }
    }
    string <- paste(string, "      </tr>\n", sep="")
  }
  
  # stars note
  if (!is.null(custom.note)) {
    note <- custom.note
  } else if (is.null(stars)) {
    note <- "\n"
  } else {
    st <- sort(stars)
    if (length(unique(st)) != length(st)) {
      stop("Duplicate elements are not allowed in the stars argument.")
    }
    if (length(st) == 4) {
      note <- paste("<sup>", star.symbol, star.symbol, star.symbol, 
          "</sup>p &lt; ", st[1], ", <sup>", star.symbol, star.symbol, 
          "</sup>p &lt; ", st[2], ", <sup>", star.symbol, "</sup>p &lt; ", st[3], 
          ", <sup>", symbol, "</sup>p &lt; ", st[4], sep="")
    } else if (length(st) == 3) {
      note <- paste("<sup>", star.symbol, star.symbol, star.symbol, 
          "</sup>p &lt; ", st[1], ", <sup>", star.symbol, star.symbol, 
          "</sup>p &lt; ", st[2], ", <sup>", star.symbol, "</sup>p &lt; ", st[3], 
          sep="")
    } else if (length(st) == 2) {
      note <- paste("<sup>", star.symbol, star.symbol, "</sup>p &lt; ", st[1], 
          ", <sup>", star.symbol, "</sup>p &lt; ", st[2], sep="")
    } else if (length(st) == 1) {
      note <- paste("<sup>", star.symbol, "</sup>p &lt; ", st[1], sep="")
    } else {
      note <- ""
    }
  }
  
  string <- paste(string, "      <tr>\n", "        <td colspan=\"", 
      (1 + length(models)), "\"><span style=\"font-size:0.8em\">", note, 
      "</span></td>\n", "      </tr>\n", sep="")
  
  # write table footer
  string <- paste(string, "    </table>\n")
  
  string <- paste(string, "  </body>\n</html>\n\n", sep="")
  
  if (is.na(file)) {
    cat(string)
  } else if (!is.character(file)) {
    stop("The 'file' argument must be a character string.")
  } else {
    sink(file)
    cat(string)
    sink()
    cat(paste("The table was written to the file '", file, "'.\n", sep=""))
  }
  if (return.string == TRUE) {
    return(string)
  }
}

