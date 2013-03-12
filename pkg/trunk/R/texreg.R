# The texreg package was written by Philip Leifeld.
# Please use the forum at http://r-forge.r-project.org/projects/texreg/ 
# for bug reports, help or feature requests.


# screenreg function
screenreg <- function(l, file=NA, single.row=FALSE, stars=c(0.001, 0.01, 0.05), 
    custom.model.names=NULL, custom.coef.names=NULL, custom.gof.names=NULL, 
    custom.note=NULL, digits=2, leading.zero=TRUE, symbol=".", override.coef=0, 
    override.se=0, override.pval=0, omit.coef=NA, reorder.coef=NULL, 
    reorder.gof=NULL, return.string=FALSE, column.spacing=2, outer.rule="=", 
    inner.rule="-", ...) {
  
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
  
  m <- customnames(m, custom.coef.names) #rename coefficients
  m <- rearrangeMatrix(m) #resort matrix and conflate duplicate entries
  m <- as.data.frame(m)
  m <- omitcoef(m, omit.coef) #remove coefficient rows matching regex
  
  modnames <- modelnames(models, custom.model.names) #use (custom) model names
  
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

texreg <- function(l, file=NA, single.row=FALSE, stars=c(0.001, 0.01, 0.05), 
    custom.model.names=NULL, custom.coef.names=NULL, custom.gof.names=NULL, 
    custom.note=NULL, digits=2, leading.zero=TRUE, symbol="\\cdot", 
    override.coef=0, override.se=0, override.pval=0, omit.coef=NA, 
    reorder.coef=NULL, reorder.gof=NULL, return.string=FALSE, bold=0.00, 
    center=TRUE, caption="Statistical models", caption.above=FALSE, 
    label="table:coefficients", booktabs=FALSE, dcolumn=FALSE, sideways=FALSE, 
    use.packages=TRUE, table=TRUE, no.margin=TRUE, scriptsize=FALSE, 
    float.pos="", ...) {
  
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
  
  m <- customnames(m, custom.coef.names) #rename coefficients
  m <- rearrangeMatrix(m) #resort matrix and conflate duplicate entries
  m <- as.data.frame(m)
  m <- omitcoef(m, omit.coef) #remove coefficient rows matching regex
  
  modnames <- modelnames(models, custom.model.names) #use (custom) model names
  
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
  if (single.row == TRUE) {
    separator <- ")"
  } else {
    separator <- "."
  }
  if (no.margin == FALSE) {
    margin.arg <- ""
  } else {
    margin.arg <- "@{}"
  }
  for (i in 2:ncol(output.matrix)) {
    if (dcolumn == FALSE) {
      string <- paste0(string, "c ")
    } else {
      if (single.row == TRUE) {
        dl <- compute.width(output.matrix[,i], left=TRUE, single.row=TRUE)
        dr <- compute.width(output.matrix[,i], left=FALSE, single.row=TRUE)
      } else {
        dl <- compute.width(output.matrix[,i], left=TRUE, single.row=FALSE)
        dr <- compute.width(output.matrix[,i], left=FALSE, single.row=FALSE)
      }
      string <- paste0(string, "D{", separator, "}{", separator, "}{", 
          dl, separator, dr, "}", margin.arg, " ")
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
        "}{l}{\\scriptsize{\\textsuperscript{***}$p<", st[1], 
        "$, \n  \\textsuperscript{**}$p<", st[2], 
        "$, \n  \\textsuperscript{*}$p<", st[3], 
        "$, \n  \\textsuperscript{$", symbol, "$}$p<", st[4], "$}}\n", sep="")
    } else if (length(st) == 3) {
      note <- paste("\\multicolumn{", length(models)+1, 
        "}{l}{\\scriptsize{\\textsuperscript{***}$p<", st[1], 
        "$, \n  \\textsuperscript{**}$p<", st[2], 
        "$, \n  \\textsuperscript{*}$p<", st[3], "$}}\n", sep="")
    } else if (length(st) == 2) {
      note <- paste("\\multicolumn{", length(models)+1, 
        "}{l}{\\scriptsize{\\textsuperscript{**}$p<", st[1], 
        "$, \n  \\textsuperscript{*}$p<", st[2], "$}}\n", sep="")
    } else if (length(st) == 1) {
      note <- paste("\\multicolumn{", length(models)+1, 
        "}{l}{\\scriptsize{\\textsuperscript{*}$p<", st[1], "$}}\n", sep="")
    } else {
      note <- ""
    }
  }
  string <- paste0(string, note, "\\end{tabular}\n")
  
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
htmlreg <- function(l, file=NA, single.row=FALSE, stars=c(0.001, 0.01, 0.05), 
    custom.model.names=NULL, custom.coef.names=NULL, custom.gof.names=NULL, 
    custom.note=NULL, digits=2, leading.zero=TRUE, symbol="&middot;", 
    override.coef=0, override.se=0, override.pval=0, omit.coef=NA, 
    reorder.coef=NULL, reorder.gof=NULL, return.string=FALSE, bold=0.00, 
    center=TRUE, caption="Statistical models", caption.above=FALSE, 
    star.symbol="*", inline.css=TRUE, doctype=TRUE, html.tag=FALSE, 
    head.tag=FALSE, body.tag=FALSE, ...) {
  
  stars <- check.stars(stars)
  
  models <- get.data(l, ...) #extract relevant coefficients, SEs, GOFs, etc.
  
  # inline CSS definitions
  if (inline.css == TRUE) {
    css.table <- " style=\"border: none;\""
    css.th <- paste0(" style=\"text-align: left; border-top: 2px solid ", 
        "black; border-bottom: 1px solid black; padding-right: 12px;\"")
    css.midrule <- " style=\"border-top: 1px solid black;\""
    css.bottomrule <- " style=\"border-bottom: 2px solid black;\""
    css.td <- " style=\"padding-right: 12px; border: none;\""
    css.caption <- " float: left;"
    css.sup <- " style=\"vertical-align: 4px;\""
  } else {
    css.table <- ""
    css.th <- ""
    css.midrule <- ""
    css.bottomrule <- ""
    css.td <- ""
    css.caption <- ""
    css.sup <- ""
  }
  
  models <- override(models, override.coef, override.se, override.pval)
  models <- tex.replace(models, type="html", style=css.sup) #TeX --> HTML
  gof.names <- get.gof(models) #extract names of GOFs
  
  # arrange coefficients and GOFs nicely in a matrix
  gofs <- aggregate.matrix(models, gof.names, custom.gof.names, digits, 
      returnobject="gofs")
  m <- aggregate.matrix(models, gof.names, custom.gof.names, digits, 
      returnobject="m")
  decimal.matrix <- aggregate.matrix(models, gof.names, custom.gof.names, 
      digits, returnobject="decimal.matrix")
  
  m <- customnames(m, custom.coef.names) #rename coefficients
  m <- rearrangeMatrix(m) #resort matrix and conflate duplicate entries
  m <- as.data.frame(m)
  m <- omitcoef(m, omit.coef) #remove coefficient rows matching regex
  
  modnames <- modelnames(models, custom.model.names) #use (custom) model names
  
  # reorder GOF and coef matrix
  m <- reorder(m, reorder.coef)
  gofs <- reorder(gofs, reorder.gof)
  decimal.matrix <- reorder(decimal.matrix, reorder.gof)
  
  # create output table with significance stars etc.
  output.matrix <- outputmatrix(m, single.row, neginfstring="-Inf", 
      leading.zero, digits, se.prefix=" (", se.suffix=")", 
      star.char=star.symbol, star.prefix=paste0("<sup", css.sup, ">"), 
      star.suffix="</sup>", stars, dcolumn=TRUE, symbol, bold, 
      bold.prefix="<b>", bold.suffix="</b>")
  
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
    doct <- paste0("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 ", 
        "Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n")
  } else {
    doct <- ""
  }
  
  # determine indentation for table
  if (html.tag == TRUE) {
    h.ind <- "  "
  } else {
    h.ind <- ""
  }
  if (body.tag == TRUE) {
    b.ind <- "  "
  } else {
    b.ind <- ""
  }
  if (head.tag == TRUE) {
    d.ind <- "  "
  } else {
    d.ind <- ""
  }
  ind <- "  "
  
  # horizontal table alignment
  if (center == FALSE) {
    tabdef <- paste0(h.ind, b.ind, "<table cellspacing=\"0\"", css.table, ">\n")
  } else {
    tabdef <- paste0(h.ind, b.ind, 
        "<table cellspacing=\"0\" align=\"center\"", css.table, ">\n")
  }
  
  # set caption
  if (is.null(caption) || !is.character(caption)) {
    stop("The caption must be provided as a (possibly empty) character vector.")
  } else if (caption != "" && caption.above == FALSE) {
    cap <- paste(h.ind, b.ind, ind, 
        "<caption align=\"bottom\" style=\"margin-top:0.3em;", css.caption, 
        "\">", caption, "</caption>\n", sep="")
  } else if (caption != "" && caption.above == TRUE) {
    cap <- paste(h.ind, b.ind, ind, 
        "<caption align=\"top\" style=\"margin-bottom:0.3em;", css.caption, 
        "\">", caption, "</caption>\n", sep="")
  } else {
    cap <- ""
  }
  
  # HTML header with CSS definitions
  string <- paste0("\n", doct)
  if (html.tag == TRUE) {
    string <- paste0(string, "<html>\n")
  }
  
  if (inline.css == TRUE) {
    css.header <- ""
  } else {
    css.header <- paste0(
        h.ind, d.ind, "<style type=\"text/css\">\n", 
        h.ind, d.ind, ind, "table {\n", 
        h.ind, d.ind, ind, ind, "border: none;\n", 
        h.ind, d.ind, ind, "}\n", 
        h.ind, d.ind, ind, "th {\n", 
        h.ind, d.ind, ind, ind, "text-align: left;\n", 
        h.ind, d.ind, ind, ind, "border-top: 2px solid black;\n", 
        h.ind, d.ind, ind, ind, "border-bottom: 1px solid black;\n", 
        h.ind, d.ind, ind, ind, "padding-right: 12px;\n", 
        h.ind, d.ind, ind, "}\n", 
        h.ind, d.ind, ind, ".midRule {\n", 
        h.ind, d.ind, ind, ind, "border-top: 1px solid black;\n", 
        h.ind, d.ind, ind, "}\n", 
        h.ind, d.ind, ind, ".bottomRule {\n", 
        h.ind, d.ind, ind, ind, "border-bottom: 2px solid black;\n", 
        h.ind, d.ind, ind, "}\n", 
        h.ind, d.ind, ind, "td {\n", 
        h.ind, d.ind, ind, ind, "padding-right: 12px;\n", 
        h.ind, d.ind, ind, ind, "border: none;\n", 
        h.ind, d.ind, ind, "}\n", 
        h.ind, d.ind, ind, "caption span {\n", 
        h.ind, d.ind, ind, ind, "float: left;\n", 
        h.ind, d.ind, ind, "}\n", 
        h.ind, d.ind, ind, "sup {\n", 
        h.ind, d.ind, ind, ind, "vertical-align: 4px;\n", 
        h.ind, d.ind, ind, "}\n", 
        h.ind, d.ind, "</style>\n"
    )
  }

  if (head.tag == TRUE) {
    string <- paste0(string, 
        h.ind, "<head>\n", 
        h.ind, d.ind, "<title>", caption, "</title>\n", 
        css.header, 
        h.ind, "</head>\n\n")
  }
  if (body.tag == TRUE) {
    string <- paste0(string, h.ind, "<body>\n")
  }
  string <- paste(
      string, 
      tabdef, 
      cap, 
      h.ind, b.ind, ind, "<tr>\n", 
      h.ind, b.ind, ind, ind, "<th", css.th, "></th>\n", 
      sep="")
  
  # specify model names (header row)
  for (i in 1:length(models)) {
    string <- paste(string, 
        h.ind, b.ind, ind, ind, "<th", css.th, "><b>", 
        modnames[i], "</b></th>\n", sep="")
  }
  string <- paste(string, h.ind, b.ind, ind, "</tr>\n", sep="")
  
  # write coefficients to string object
  for (i in 1:(length(output.matrix[,1])-length(gof.names))) {
    string <- paste(string, h.ind, b.ind, ind, "<tr>\n", sep="")
    for (j in 1:length(output.matrix[1,])) {
      string <- paste0(string, h.ind, b.ind, ind, ind, "<td", css.td, ">", 
      output.matrix[i,j], "</td>\n")
    }
    string <- paste(string, h.ind, b.ind, ind, "</tr>\n", sep="")
  }
  
  for (i in (length(output.matrix[,1])-(length(gof.names)-1)):
      (length(output.matrix[,1]))) {
    string <- paste(string, h.ind, b.ind, ind, "<tr>\n", sep="")
    for (j in 1:length(output.matrix[1,])) {
      if (i == length(output.matrix[,1])-(length(gof.names)-1)) {
        if (inline.css == TRUE) {
          mr <- css.midrule
        } else {
          mr <- " class=\"midRule\""
        }
        string <- paste(string, h.ind, b.ind, ind, ind, 
            "<td", mr, ">", output.matrix[i,j], "</td>\n", sep="")
      } else if (i == length(output.matrix[,1])) {
        if (inline.css == TRUE) {
          br <- css.bottomrule
        } else {
          br <- " class=\"bottomRule\""
        }
        string <- paste(string, h.ind, b.ind, ind, ind, 
            "<td", br, ">", output.matrix[i,j], "</td>\n", sep="")
      } else {
        string <- paste(string, h.ind, b.ind, ind, ind, "<td", css.td, ">", 
            output.matrix[i,j], "</td>\n", sep="")
      }
    }
    string <- paste(string, h.ind, b.ind, ind, "</tr>\n", sep="")
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
      note <- paste("<sup", css.sup, ">", star.symbol, star.symbol, star.symbol,
          "</sup>p &lt; ", st[1], ", <sup", css.sup, ">", 
          star.symbol, star.symbol, "</sup", css.sup, ">p &lt; ", st[2], 
          ", <sup", css.sup, ">", star.symbol, "</sup>p &lt; ", 
          st[3], ", <sup", css.sup, ">", symbol, "</sup>p &lt; ", 
          st[4], sep="")
    } else if (length(st) == 3) {
      note <- paste("<sup", css.sup, ">", star.symbol, star.symbol, star.symbol,
          "</sup>p &lt; ", st[1], ", <sup", css.sup, ">", star.symbol, 
          star.symbol, "</sup>p &lt; ", st[2], ", <sup", css.sup, ">", 
          star.symbol, "</sup>p &lt; ", st[3], sep="")
    } else if (length(st) == 2) {
      note <- paste("<sup", css.sup, ">", star.symbol, star.symbol, 
          "</sup>p &lt; ", st[1], ", <sup", css.sup, ">", star.symbol, 
          "</sup>p &lt; ", st[2], sep="")
    } else if (length(st) == 1) {
      note <- paste0("<sup", css.sup, ">", star.symbol, "</sup>p &lt; ", st[1])
    } else {
      note <- ""
    }
  }
  
  string <- paste(string, h.ind, b.ind, ind, "<tr>\n", h.ind, b.ind, ind, ind, 
      "<td", css.td, " colspan=\"", (1 + length(models)), 
      "\"><span style=\"font-size:0.8em\">", note, "</span></td>\n", h.ind, 
      b.ind, ind, "</tr>\n", sep="")
  
  # write table footer
  string <- paste0(string, h.ind, b.ind, "</table>\n")
  if (body.tag == TRUE) {
    string <- paste0(string, h.ind, "</body>\n")
  }
  if (html.tag == TRUE) {
    string <- paste0(string, "</html>\n\n")
  } else {
    string <- paste0(string, "\n")
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

