# The texreg package was written by Philip Leifeld.
# Please use the forum at http://r-forge.r-project.org/projects/texreg/ 
# for bug reports, help or feature requests.


# screenreg function
screenreg <- function(l, file = NA, single.row = FALSE, 
    stars = c(0.001, 0.01, 0.05), custom.model.names = NULL, 
    custom.coef.names = NULL, custom.gof.names = NULL, custom.note = NULL, 
    digits = 2, leading.zero = TRUE, symbol = ".", override.coef = 0, 
    override.se = 0, override.pval = 0, omit.coef = NA, reorder.coef = NULL, 
    reorder.gof = NULL, return.string = FALSE, ci.force = FALSE,
    ci.force.level = 0.95, ci.test = 0, column.spacing = 2, outer.rule = "=", 
    inner.rule = "-", ...) {
  
  stars <- check.stars(stars)
  
  models <- get.data(l, ...)  #extract relevant coefficients, SEs, GOFs, etc.
  models <- override(models, override.coef, override.se, override.pval)
  models <- tex.replace(models, type = "screen")  #convert TeX code to text code
  models <- ciforce(models, ci.force = ci.force, ci.level = ci.force.level)
  gof.names <- get.gof(models)  #extract names of GOFs
  
  # arrange coefficients and GOFs nicely in a matrix
  gofs <- aggregate.matrix(models, gof.names, custom.gof.names, digits, 
      returnobject = "gofs")
  m <- aggregate.matrix(models, gof.names, custom.gof.names, digits, 
      returnobject = "m")
  decimal.matrix <- aggregate.matrix(models, gof.names, custom.gof.names, 
      digits, returnobject = "decimal.matrix")
  
  m <- customnames(m, custom.coef.names)  #rename coefficients
  m <- rearrangeMatrix(m)  #resort matrix and conflate duplicate entries
  m <- as.data.frame(m)
  m <- omitcoef(m, omit.coef)  #remove coefficient rows matching regex
  
  modnames <- modelnames(models, custom.model.names)  #use (custom) model names
  
  # reorder GOF and coef matrix
  m <- reorder(m, reorder.coef)
  gofs <- reorder(gofs, reorder.gof)
  decimal.matrix <- reorder(decimal.matrix, reorder.gof)
  
  # create output table with significance stars etc.
  ci <- logical()
  for (i in 1:length(models)) {
    if (length(models[[i]]@se) == 0) {
      ci[i] <- TRUE
    } else {
      ci[i] <- FALSE
    }
  }
  output.matrix <- outputmatrix(m, single.row, neginfstring = "-Inf", 
      leading.zero, digits, se.prefix = " (", se.suffix = ")", 
      star.prefix = " ", star.suffix = "", star.char = "*", stars, 
      dcolumn = TRUE, symbol = symbol, bold = 0, bold.prefix = "", 
      bold.suffix = "", ci = ci, ci.test = ci.test)
  
  # create GOF matrix (the lower part of the final output matrix)
  gof.matrix <- gofmatrix(gofs, decimal.matrix, dcolumn = TRUE, leading.zero, 
      digits)
  
  # combine the coefficient and gof matrices vertically
  output.matrix <- rbind(output.matrix, gof.matrix)
  
  # reformat output matrix and add spaces
  if (ncol(output.matrix) == 2) {
    temp <- matrix(format.column(output.matrix[, -1], single.row = single.row, 
        digits = digits))
  } else {
    temp <- apply(output.matrix[, -1], 2, format.column, 
        single.row = single.row, digits = digits)
  }
  output.matrix <- cbind(output.matrix[, 1], temp)
  output.matrix <- rbind(c("", modnames), output.matrix)
  for (i in 1:ncol(output.matrix)) {
    output.matrix[, i] <- fill.spaces(output.matrix[, i])
  }
  
  string <- "\n"
  
  # horizontal rule above the table
  table.width <- sum(nchar(output.matrix[1, ])) + 
      (ncol(output.matrix) - 1) * column.spacing
  if (class(outer.rule) != "character") {
    stop("outer.rule must be a character.")
  } else if (nchar(outer.rule) > 1) {
    stop("outer.rule must be a character of maximum length 1.")
  } else if (outer.rule == "") {
    o.rule <- ""
  } else {
    o.rule <- paste(rep(outer.rule, table.width), collapse = "")
    string <- paste0(string, o.rule, "\n")
  }
  
  # specify model names
  spacing <- paste(rep(" ", column.spacing), collapse = "")
  string <- paste(string, output.matrix[1, 1], sep = "")
  for (i in 2:ncol(output.matrix)) {
    string <- paste0(string, spacing, output.matrix[1, i])
  }
  string <- paste0(string, "\n")
  
  # mid rule 1
  if (class(inner.rule) != "character") {
    stop("inner.rule must be a character.")
  } else if (nchar(inner.rule) > 1) {
    stop("inner.rule must be a character of maximum length 1.")
  } else if (inner.rule == "") {
    i.rule <- ""
  } else {
    i.rule <- paste(rep(inner.rule, table.width), collapse = "")
    string <- paste0(string, i.rule, "\n")
  }
  
  # write coefficients
  for (i in 2:(length(output.matrix[, 1]) - length(gof.names))) {
    for (j in 1:length(output.matrix[1, ])) {
      string <- paste0(string, output.matrix[i,j])
      if (j == length(output.matrix[1, ])) {
        string <- paste0(string, "\n")
      } else {
        string <- paste0(string, spacing)
      }
    }
  }
  
  # mid rule 2
  if (inner.rule != "") {
    string <- paste0(string, i.rule, "\n")
  }
  
  # write GOF part of the output matrix
  for (i in (length(output.matrix[, 1]) - (length(gof.names) - 1)):
      (length(output.matrix[, 1]))) {
    for (j in 1:length(output.matrix[1, ])) {
      string <- paste0(string, output.matrix[i, j])
      if (j == length(output.matrix[1, ])) {
        string <- paste0(string, "\n")
      } else {
        string <- paste0(string, spacing)
      }
    }
  }
  
  # write table footer
  if (outer.rule != "") {
    string <- paste0(string, o.rule, "\n")
  }
  
  # stars note
  if (is.null(stars)) {
    snote <- ""
  } else if (any(ci == FALSE)) {
    st <- sort(stars)
    if (length(unique(st)) != length(st)) {
      stop("Duplicate elements are not allowed in the stars argument.")
    }
    if (length(st) == 4) {
      snote <- paste0("*** p < ", st[1], ", ** p < ", st[2], ", * p < ", st[3], 
          ", ", symbol, " p < ", st[4])
    } else if (length(st) == 3) {
      snote <- paste0("*** p < ", st[1], ", ** p < ", st[2], ", * p < ", st[3])
    } else if (length(st) == 2) {
      snote <- paste0("** p < ", st[1], ", * p < ", st[2])
    } else if (length(st) == 1) {
      snote <- paste0("* p < ", st)
    } else {
      snote <- ""
    }
    if (is.numeric(ci.test) && !is.na(ci.test) && nchar(snote) > 0 && any(ci)) {
      snote <- paste(snote, "(or", ci.test, "outside the confidence interval).")
    } else if (is.numeric(ci.test) && !is.na(ci.test) && any(ci)) {
      snote <- paste("*", ci.test, "outside the confidence interval")
    }
  } else if (is.numeric(ci.test) && !is.na(ci.test)) {
    snote <- paste("*", ci.test, "outside the confidence interval")
  } else {
    snote <- ""
  }
  if (is.null(custom.note)) {
    note <- paste0(snote, "\n\n")
  } else if (custom.note == "") {
    note <- "\n"
  } else {
    note <- paste0(custom.note, "\n\n")
    note <- gsub("%stars", snote, note)
  }
  string <- paste0(string, note)
  
  #write to file
  if (is.na(file)) {
    cat(string)
  } else if (!is.character(file)) {
    stop("The 'file' argument must be a character string.")
  } else {
    sink(file)
    cat(string)
    sink()
    cat(paste0("The table was written to the file '", file, "'.\n"))
  }
  
  if (return.string == TRUE) {
    return(string)
  }
}


# texreg function

texreg <- function(l, file = NA, single.row = FALSE, 
    stars = c(0.001, 0.01, 0.05), custom.model.names = NULL, 
    custom.coef.names = NULL, custom.gof.names = NULL, custom.note = NULL, 
    digits = 2, leading.zero = TRUE, symbol = "\\cdot", override.coef = 0, 
    override.se = 0, override.pval = 0, omit.coef = NA, reorder.coef = NULL, 
    reorder.gof = NULL, return.string = FALSE, ci.force = FALSE,
    ci.force.level = 0.95, ci.test = 0, bold = 0.00, center = TRUE, 
    caption = "Statistical models", caption.above = FALSE, 
    label = "table:coefficients", booktabs = FALSE, dcolumn = FALSE, 
    sideways = FALSE, use.packages = TRUE, table = TRUE, no.margin = TRUE, 
    scriptsize = FALSE, float.pos = "", ...) {
  
  stars <- check.stars(stars)
  
  #check dcolumn vs. bold
  if (dcolumn == TRUE && bold > 0) {
    dcolumn <- FALSE
    msg <- paste("The dcolumn package and the bold argument cannot be used at", 
        "the same time. Switching off dcolumn.")
    if (stars == TRUE) {
      warning(paste(msg, "You should also consider setting stars = FALSE."))
    } else {
      warning(msg)
    }
  }
  
  models <- get.data(l, ...)  #extract relevant coefficients, SEs, GOFs, etc.
  gof.names <- get.gof(models)  #extract names of GOFs
  models <- override(models, override.coef, override.se, override.pval)
  models <- ciforce(models, ci.force = ci.force, ci.level = ci.force.level)
  
  # arrange coefficients and GOFs nicely in a matrix
  gofs <- aggregate.matrix(models, gof.names, custom.gof.names, digits, 
      returnobject = "gofs")
  m <- aggregate.matrix(models, gof.names, custom.gof.names, digits, 
      returnobject = "m")
  decimal.matrix <- aggregate.matrix(models, gof.names, custom.gof.names, 
      digits, returnobject = "decimal.matrix")
  
  m <- customnames(m, custom.coef.names)  #rename coefficients
  m <- rearrangeMatrix(m)  #resort matrix and conflate duplicate entries
  m <- as.data.frame(m)
  m <- omitcoef(m, omit.coef)  #remove coefficient rows matching regex
  
  modnames <- modelnames(models, custom.model.names)  #use (custom) model names
  
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
  ci <- logical()
  for (i in 1:length(models)) {
    if (length(models[[i]]@se) == 0) {
      ci[i] <- TRUE
    } else {
      ci[i] <- FALSE
    }
  }
  output.matrix <- outputmatrix(m, single.row, 
      neginfstring = "\\multicolumn{1}{c}{$-$Inf}", leading.zero, digits, 
      se.prefix = " \\; (", se.suffix = ")", star.prefix = "^{", 
      star.suffix = "}", star.char = "*", stars, dcolumn = dcolumn, 
      symbol, bold, bold.prefix = "\\textbf{", bold.suffix = "}", ci = ci, 
      semicolon = ";\\ ", ci.test = ci.test)
  
  # create GOF matrix (the lower part of the final output matrix)
  gof.matrix <- gofmatrix(gofs, decimal.matrix, dcolumn = TRUE, leading.zero, 
      digits)
  
  # combine the coefficient and gof matrices vertically
  output.matrix <- rbind(output.matrix, gof.matrix)
  
  string <- ""
  
  # write table header
  string <- paste0(string, "\n")
  if (use.packages == TRUE) {
    if (sideways == TRUE & table == TRUE) {
      string <- paste0(string, "\\usepackage{rotating}\n")
    }
    if (booktabs == TRUE) {
      string <- paste0(string, "\\usepackage{booktabs}\n")
    }
    if (dcolumn == TRUE) {
      string <- paste0(string, "\\usepackage{dcolumn}\n")
    }
    if (dcolumn == TRUE || booktabs == TRUE || sideways == TRUE) {
      cat("\n")
    }
  }
  if (table == TRUE) {
    if (sideways == TRUE) {
      t <- "sideways"
    } else {
      t <- ""
    }
    if ( float.pos == "") {
      string <- paste0(string, "\\begin{", t, "table}\n")
    } else {
      string <- paste0(string, "\\begin{", t, "table}[", float.pos, "]\n")
    }
    if (caption.above == TRUE) {
      string <- paste0(string, "\\caption{", caption, "}\n")
    }
    if (center == TRUE) {
      string <- paste0(string, "\\begin{center}\n")
    }
    if (scriptsize == TRUE) {
      string <- paste0(string, "\\scriptsize\n")
    }
  }
  string <- paste0(string, "\\begin{tabular}{l ")
  
  #define columns of the table
  if (no.margin == FALSE) {
    margin.arg <- ""
  } else {
    margin.arg <- "@{}"
  }
  for (i in 2:ncol(output.matrix)) {
    if (single.row == TRUE) {
      if (ci[i - 1] == FALSE) {
        separator <- ")"
      } else {
        separator <- "]"
      }
    } else {
      separator <- "."
    }
    if (dcolumn == FALSE) {
      string <- paste0(string, "c ")
    } else {
      if (single.row == TRUE) {
        dl <- compute.width(output.matrix[, i], left = TRUE, single.row = TRUE, 
            bracket = separator)
        dr <- compute.width(output.matrix[, i], left = FALSE, single.row = TRUE,
            bracket = separator)
      } else {
        dl <- compute.width(output.matrix[, i], left = TRUE, single.row = FALSE,
            bracket = separator)
        dr <- compute.width(output.matrix[, i], left = FALSE, 
            single.row = FALSE, bracket = separator)
      }
      string <- paste0(string, "D{", separator, "}{", separator, "}{", 
          dl, separator, dr, "}", margin.arg, " ")
    }
  }
  
  # horizontal rule above the table
  if (booktabs == TRUE) {
    string <- paste0(string, "}\n", "\\toprule\n")
  } else {
    string <- paste0(string, "}\n", "\\hline\n")
  }
  
  # specify model names
  for (k in 1:lab.length) {
    string <- paste0(string, " ")
  }

  if (dcolumn == TRUE) {
    for (i in 1:length(models)) {
      string <- paste0(string, " & \\multicolumn{1}{c}{", modnames[i], "}")
    }
  } else {
    for (i in 1:length(models)) {
      string <- paste0(string, " & ", modnames[i])
    }
  }
  
  # horizontal rule between coefficients and goodness-of-fit block
  if (booktabs == TRUE) {
    string <- paste0(string, " \\\\\n", "\\midrule\n")
  } else {
    string <- paste0(string, " \\\\\n", "\\hline\n")
  }
  
  # fill with spaces
  max.lengths <- numeric(length(output.matrix[1, ]))
  for (i in 1:length(output.matrix[1, ])) {
    max.length <- 0
    for (j in 1:length(output.matrix[, 1])) {
      if (nchar(output.matrix[j, i]) > max.length) {
        max.length <- nchar(output.matrix[j, i])
      }
    }
    max.lengths[i] <- max.length
  }
  for (i in 1:length(output.matrix[, 1])) {
    for (j in 1:length(output.matrix[1, ])) {
      nzero <- max.lengths[j] - nchar(output.matrix[i, j])
      zeros <- rep(" ", nzero)
      zeros <- paste(zeros, collapse = "")
      output.matrix[i, j] <- paste0(output.matrix[i, j], zeros)
    }
  }
  
  # write coefficients to string object
  for (i in 1:(length(output.matrix[, 1]) - length(gof.names))) {
    for (j in 1:length(output.matrix[1, ])) {
      string <- paste0(string, output.matrix[i, j])
      if (j == length(output.matrix[1, ])) {
        string <- paste0(string, " \\\\\n")
      } else {
        string <- paste0(string, " & ")
      }
    }
  }
  
  if (booktabs == TRUE) {
    string <- paste0(string, "\\midrule\n")
  } else {
    string <- paste0(string, "\\hline\n")
  }
  
  for (i in (length(output.matrix[, 1]) - (length(gof.names) - 1)):
      (length(output.matrix[, 1]))) {
    for (j in 1:length(output.matrix[1, ])) {
      string <- paste0(string, output.matrix[i, j])
      if (j == length(output.matrix[1, ])) {
        string <- paste0(string, " \\\\\n")
      } else {
        string <- paste0(string, " & ")
      }
    }
  }
  
  # write table footer
  if (booktabs == TRUE) {
    string <- paste0(string, "\\bottomrule\n")
  } else {
    string <- paste0(string, "\\hline\n")
  }
  
  # stars note
  if (is.null(stars)) {
    snote <- ""
  } else if (any(ci == FALSE)) {
    st <- sort(stars)
    if (length(unique(st)) != length(st)) {
      stop("Duplicate elements are not allowed in the stars argument.")
    }
    if (length(st) == 4) {
      snote <- paste0("\\textsuperscript{***}$p<", st[1], 
        "$, \n  \\textsuperscript{**}$p<", st[2], 
        "$, \n  \\textsuperscript{*}$p<", st[3], 
        "$, \n  \\textsuperscript{$", symbol, "$}$p<", st[4], "$")
    } else if (length(st) == 3) {
      snote <- paste0("\\textsuperscript{***}$p<", st[1], 
        "$, \n  \\textsuperscript{**}$p<", st[2], 
        "$, \n  \\textsuperscript{*}$p<", st[3], "$")
    } else if (length(st) == 2) {
      snote <- paste0("\\textsuperscript{**}$p<", st[1], 
        "$, \n  \\textsuperscript{*}$p<", st[2], "$")
    } else if (length(st) == 1) {
      snote <- paste0("\\textsuperscript{*}$p<", st[1], "$")
    } else {
      snote <- ""
    }
    if (is.numeric(ci.test) && !is.na(ci.test) && nchar(snote) > 0 && any(ci)) {
      snote <- paste(snote, "(or", ci.test, "outside the confidence interval).")
    } else if (is.numeric(ci.test) && !is.na(ci.test) && any(ci)) {
      snote <- paste("\\textsuperscript{*}", ci.test, 
          "outside the confidence interval")
    }
  } else if (is.numeric(ci.test) && !is.na(ci.test)) {
    snote <- paste("\\textsuperscript{*}", ci.test, 
        "outside the confidence interval")
  } else {
    snote <- ""
  }
  if (is.null(custom.note)) {
    note <- paste0("\\multicolumn{", length(models) + 1, 
        "}{l}{\\scriptsize{", snote, "}}\n")
  } else if (custom.note == "") {
    note <- ""
  } else {
    note <- paste0("\\multicolumn{", length(models) + 1, 
        "}{l}{\\scriptsize{", custom.note, "}}\n")
    note <- gsub("%stars", snote, note)
  }
  string <- paste0(string, note, "\\end{tabular}\n")
  
  if (table == TRUE) {
    if (scriptsize == TRUE) {
      string <- paste0(string, "\\normalsize\n")
    }
    if (caption.above == FALSE) {
      string <- paste0(string, "\\caption{", caption, "}\n")
    }
    string <- paste0(string, "\\label{", label, "}\n")
    if (center == TRUE) {
      string <- paste0(string, "\\end{center}\n")
    }
    if (sideways == TRUE) {
      t <- "sideways"
    } else {
      t <- ""
    }
    string <- paste0(string, "\\end{", t, "table}\n\n")
  }
  
  if (is.na(file)) {
    cat(string)
  } else if (!is.character(file)) {
    stop("The 'file' argument must be a character string.")
  } else {
    sink(file)
    cat(string)
    sink()
    cat(paste0("The table was written to the file '", file, "'.\n"))
  }
  if (return.string == TRUE) {
    return(string)
  }
}


# htmlreg function
htmlreg <- function(l, file = NA, single.row = FALSE, 
    stars = c(0.001, 0.01, 0.05), custom.model.names = NULL, 
    custom.coef.names = NULL, custom.gof.names = NULL, custom.note = NULL, 
    digits = 2, leading.zero = TRUE, symbol = "&middot;", override.coef = 0, 
    override.se = 0, override.pval = 0, omit.coef = NA, reorder.coef = NULL, 
    reorder.gof = NULL, return.string = FALSE, ci.force = FALSE,
    ci.force.level = 0.95, ci.test = 0, bold = 0.00, center = TRUE, 
    caption = "Statistical models", caption.above = FALSE, star.symbol = "*", 
    inline.css = TRUE, doctype = TRUE, html.tag = FALSE, head.tag = FALSE, 
    body.tag = FALSE, ...) {
  
  stars <- check.stars(stars)
  
  models <- get.data(l, ...)  #extract relevant coefficients, SEs, GOFs, etc.
  
  # inline CSS definitions
  if (inline.css == TRUE) {
    css.table <- " style=\"border: none;\""
    css.th <- paste0(" style=\"text-align: left; border-top: 2px solid ", 
        "black; border-bottom: 1px solid black; padding-right: 12px;\"")
    css.midrule <- " style=\"border-top: 1px solid black;\""
    css.bottomrule <- " style=\"border-bottom: 2px solid black;\""
    css.td <- " style=\"padding-right: 12px; border: none;\""
    css.caption <- ""
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
  models <- tex.replace(models, type = "html", style = css.sup)  # TeX --> HTML
  models <- ciforce(models, ci.force = ci.force, ci.level = ci.force.level)
  gof.names <- get.gof(models)  # extract names of GOFs
  
  # arrange coefficients and GOFs nicely in a matrix
  gofs <- aggregate.matrix(models, gof.names, custom.gof.names, digits, 
      returnobject = "gofs")
  m <- aggregate.matrix(models, gof.names, custom.gof.names, digits, 
      returnobject = "m")
  decimal.matrix <- aggregate.matrix(models, gof.names, custom.gof.names, 
      digits, returnobject = "decimal.matrix")
  
  m <- customnames(m, custom.coef.names)  # rename coefficients
  m <- rearrangeMatrix(m)  # resort matrix and conflate duplicate entries
  m <- as.data.frame(m)
  m <- omitcoef(m, omit.coef)  # remove coefficient rows matching regex
  
  modnames <- modelnames(models, custom.model.names)  # use (custom) model names
  
  # reorder GOF and coef matrix
  m <- reorder(m, reorder.coef)
  gofs <- reorder(gofs, reorder.gof)
  decimal.matrix <- reorder(decimal.matrix, reorder.gof)
  
  # create output table with significance stars etc.
  ci <- logical()
  for (i in 1:length(models)) {
    if (length(models[[i]]@se) == 0) {
      ci[i] <- TRUE
    } else {
      ci[i] <- FALSE
    }
  }
  output.matrix <- outputmatrix(m, single.row, neginfstring = "-Inf", 
      leading.zero, digits, se.prefix = " (", se.suffix = ")", 
      star.char = star.symbol, star.prefix = paste0("<sup", css.sup, ">"), 
      star.suffix = "</sup>", stars, dcolumn = TRUE, symbol, bold = bold, 
      bold.prefix = "<b>", bold.suffix = "</b>", ci = ci, ci.test = ci.test)
  
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
    cap <- paste0(h.ind, b.ind, ind, 
        "<caption align=\"bottom\" style=\"margin-top:0.3em;", css.caption, 
        "\">", caption, "</caption>\n")
  } else if (caption != "" && caption.above == TRUE) {
    cap <- paste0(h.ind, b.ind, ind, 
        "<caption align=\"top\" style=\"margin-bottom:0.3em;", css.caption, 
        "\">", caption, "</caption>\n")
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
  string <- paste0(
      string, 
      tabdef, 
      cap, 
      h.ind, b.ind, ind, "<tr>\n", 
      h.ind, b.ind, ind, ind, "<th", css.th, "></th>\n"
  )
  
  # specify model names (header row)
  for (i in 1:length(models)) {
    string <- paste0(string, 
        h.ind, b.ind, ind, ind, "<th", css.th, "><b>", modnames[i], 
        "</b></th>\n")
  }
  string <- paste0(string, h.ind, b.ind, ind, "</tr>\n")
  
  # write coefficients to string object
  for (i in 1:(length(output.matrix[, 1]) - length(gof.names))) {
    string <- paste0(string, h.ind, b.ind, ind, "<tr>\n")
    for (j in 1:length(output.matrix[1, ])) {
      string <- paste0(string, h.ind, b.ind, ind, ind, "<td", css.td, ">", 
      output.matrix[i,j], "</td>\n")
    }
    string <- paste0(string, h.ind, b.ind, ind, "</tr>\n")
  }
  
  for (i in (length(output.matrix[, 1]) - (length(gof.names) - 1)):
      (length(output.matrix[, 1]))) {
    string <- paste0(string, h.ind, b.ind, ind, "<tr>\n")
    for (j in 1:length(output.matrix[1, ])) {
      if (i == length(output.matrix[, 1]) - (length(gof.names) - 1)) {
        if (inline.css == TRUE) {
          mr <- css.midrule
        } else {
          mr <- " class=\"midRule\""
        }
        string <- paste0(string, h.ind, b.ind, ind, ind, 
            "<td", mr, ">", output.matrix[i,j], "</td>\n")
      } else if (i == length(output.matrix[, 1])) {
        if (inline.css == TRUE) {
          br <- css.bottomrule
        } else {
          br <- " class=\"bottomRule\""
        }
        string <- paste0(string, h.ind, b.ind, ind, ind, 
            "<td", br, ">", output.matrix[i,j], "</td>\n")
      } else {
        string <- paste0(string, h.ind, b.ind, ind, ind, "<td", css.td, ">", 
            output.matrix[i,j], "</td>\n")
      }
    }
    string <- paste0(string, h.ind, b.ind, ind, "</tr>\n")
  }
  
  # stars note
  if (is.null(stars)) {
    snote <- ""
  } else if (any(ci == FALSE)) {
    st <- sort(stars)
    if (length(unique(st)) != length(st)) {
      stop("Duplicate elements are not allowed in the stars argument.")
    }
    if (length(st) == 4) {
      snote <- paste0("<sup", css.sup, ">", star.symbol, star.symbol, 
          star.symbol, "</sup>p &lt; ", st[1], ", <sup", css.sup, ">", 
          star.symbol, star.symbol, "</sup", css.sup, ">p &lt; ", st[2], 
          ", <sup", css.sup, ">", star.symbol, "</sup>p &lt; ", 
          st[3], ", <sup", css.sup, ">", symbol, "</sup>p &lt; ", st[4])
    } else if (length(st) == 3) {
      snote <- paste0("<sup", css.sup, ">", star.symbol, star.symbol, 
          star.symbol, "</sup>p &lt; ", st[1], ", <sup", css.sup, ">", 
          star.symbol, star.symbol, "</sup>p &lt; ", st[2], ", <sup", css.sup, 
          ">", star.symbol, "</sup>p &lt; ", st[3])
    } else if (length(st) == 2) {
      snote <- paste0("<sup", css.sup, ">", star.symbol, star.symbol, 
          "</sup>p &lt; ", st[1], ", <sup", css.sup, ">", star.symbol, 
          "</sup>p &lt; ", st[2])
    } else if (length(st) == 1) {
      snote <- paste0("<sup", css.sup, ">", star.symbol, "</sup>p &lt; ", st[1])
    } else {
      snote <- ""
    }
    if (is.numeric(ci.test) && !is.na(ci.test) && nchar(snote) > 0 && any(ci)) {
      snote <- paste(snote, "(or", ci.test, "outside the confidence interval).")
    } else if (is.numeric(ci.test) && !is.na(ci.test) && any(ci)) {
      snote <- paste0("<sup>", star.symbol, "</sup> ", ci.test, 
          " outside the confidence interval")
    }
  } else if (is.numeric(ci.test) && !is.na(ci.test)) {
    snote <- paste0("<sup>", star.symbol, "</sup> ", ci.test, 
        " outside the confidence interval")
  } else {
    snote <- ""
  }
  if (is.null(custom.note)) {
    note <- snote
  } else if (custom.note == "") {
    note <- ""
  } else {
    note <- custom.note
    note <- gsub("%stars", snote, note)
  }
  string <- paste0(string, h.ind, b.ind, ind, "<tr>\n", h.ind, b.ind, ind, ind, 
      "<td", css.td, " colspan=\"", (1 + length(models)), 
      "\"><span style=\"font-size:0.8em\">", note, "</span></td>\n", h.ind, 
      b.ind, ind, "</tr>\n")
  
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
    cat(paste0("The table was written to the file '", file, "'.\n"))
  }
  if (return.string == TRUE) {
    return(string)
  }
}

