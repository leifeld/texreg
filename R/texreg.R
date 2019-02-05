# The texreg package was written by Philip Leifeld.
# Please use the issue tracker at http://github.com/leifeld/texreg
# for bug reports, help, or feature requests.

# matrixreg function
matrixreg <- function(l,
                      single.row = FALSE,
                      stars = c(0.001, 0.01, 0.05),
                      custom.model.names = NULL,
                      custom.coef.names = NULL,
                      custom.coef.map = NULL,
                      custom.gof.names = NULL,
                      digits = 2,
                      leading.zero = TRUE,
                      star.symbol = '*',
                      symbol = ".",
                      override.coef = 0,
                      override.se = 0,
                      override.pvalues = 0,
                      override.ci.low = 0,
                      override.ci.up = 0,
                      omit.coef = NULL,
                      reorder.coef = NULL,
                      reorder.gof = NULL,
                      ci.force = "null",
                      ci.force.level = 0.95,
                      ci.test = 0,
                      bold = 0,
                      groups = NULL,
                      custom.columns = NULL,
                      custom.col.pos = NULL,
                      dcolumn = TRUE,
                      ...) {
  
  # unnamed arguments to environment
  dots <- list(...)

  # argument for internal use
  if (!'output.type' %in% names(dots)) {
    dots[['output.type']] = 'ascii'
  }

  # extract
  models <- get.data(l, ...) # extract relevant coefficients, SEs, GOFs etc.
  models <- override(models, override.coef, override.se, override.pvalues, 
      override.ci.low, override.ci.up)
  if (dots$output.type != 'latex') {
    models <- tex.replace(models, type = "screen")  # convert TeX code to text
  }
  models <- ciforce(models, ci.force = ci.force, ci.level = ci.force.level)
  gof.names <- get.gof(models)  #extract names of GOFs
  models <- correctDuplicateCoefNames(models)
  
  # arrange coefficients and GOFs nicely in a matrix
  gofs <- aggregate.matrix(models, gof.names, custom.gof.names, digits, 
      returnobject = "gofs")
  m <- aggregate.matrix(models, gof.names, custom.gof.names, digits, 
      returnobject = "m")
  decimal.matrix <- aggregate.matrix(models, gof.names, custom.gof.names, 
      digits, returnobject = "decimal.matrix")
  
  if (!is.null(custom.coef.map)) {
    m <- custommap(m, custom.coef.map)
  } else {
    m <- omit_rename(m,
                     omit.coef = omit.coef,
                     custom.coef.names = custom.coef.names)
  }
  m <- rearrangeMatrix(m)  #resort matrix and conflate duplicate entries
  m <- as.data.frame(m)

  mod.names <- modelnames(l, models, custom.model.names)  # model names
  
  # reorder GOF and coef matrix
  m <- reorder(m, reorder.coef)
  gofs <- reorder(gofs, reorder.gof)
  decimal.matrix <- reorder(decimal.matrix, reorder.gof)
  
  # create output table with significance stars etc.
  ci <- logical()
  for (i in 1:length(models)) {
    if (length(models[[i]]@se) == 0 && length(models[[i]]@ci.up) > 0) {
      ci[i] <- TRUE
    } else {
      ci[i] <- FALSE
    }
  }
  
  # output matrix
  if (dots$output.type == 'ascii') {
    output.matrix <- outputmatrix(m, single.row, neginfstring = "-Inf", 
        posinfstring = "Inf", leading.zero, digits, 
        se.prefix = " (", se.suffix = ")", star.prefix = " ", star.suffix = "", 
        stars, dcolumn = TRUE, star.symbol = star.symbol, symbol = symbol, 
        bold = bold, bold.prefix = "", bold.suffix = "", ci = ci, 
        ci.test = ci.test)
  } else if (dots$output.type == 'latex') {
    output.matrix <- outputmatrix(m, single.row, 
        neginfstring = "\\multicolumn{1}{c}{$-\\infty$}", 
        posinfstring = "\\multicolumn{1}{c}{$\\infty$}", leading.zero, digits, 
        se.prefix = " \\; (", se.suffix = ")", star.prefix = "^{", 
        star.suffix = "}", stars, dcolumn = dcolumn, star.symbol = star.symbol,
        symbol = symbol, bold = bold, bold.prefix = "\\mathbf{", 
        bold.suffix = "}", ci = ci, semicolon = ";\\ ", ci.test = ci.test, rowLabelType = 'latex')
  } else if (dots$output.type == 'html') {
    output.matrix <- outputmatrix(m, single.row, neginfstring = "-Inf", 
        posinfstring = "Inf", leading.zero, digits, 
        se.prefix = " (", se.suffix = ")", star.symbol = star.symbol, 
        star.prefix = paste0("<sup", dots$css.sup, ">"), 
        star.suffix = "</sup>", stars, dcolumn = TRUE, symbol = symbol, 
        bold = bold, bold.prefix = "<b>", bold.suffix = "</b>", ci = ci, 
        ci.test = ci.test)
  }

  # grouping
  if (dots$output.type == 'latex') {
      output.matrix <- grouping(output.matrix, groups, indentation = "    ", 
                                single.row = single.row, prefix = "", suffix = "", rowLabelType = 'latex')
  } else {
      output.matrix <- grouping(output.matrix, groups, indentation = "    ", 
                                single.row = single.row, prefix = "", suffix = "", rowLabelType = 'text')
  }
  
  # create GOF matrix (the lower part of the final output matrix)
  if (dots$output.type == 'latex') {
      gof.matrix <- gofmatrix(gofs, decimal.matrix, dcolumn = TRUE, leading.zero, 
                              digits, rowLabelType = 'latex')
  } else {
      gof.matrix <- gofmatrix(gofs, decimal.matrix, dcolumn = TRUE, leading.zero, 
                              digits, rowLabelType = 'text')
  } 
  
  # combine the coefficient and gof matrices vertically
  coef.names <- row.names(output.matrix)
  output.matrix <- rbind(output.matrix, gof.matrix)
  
  # reformat output matrix
  if (ncol(output.matrix) == 2) {
    temp <- matrix(format.column(matrix(output.matrix[, -1]), 
        single.row = single.row, digits = digits))
  } else {
    temp <- apply(output.matrix[, -1], 2, format.column, 
        single.row = single.row, digits = digits)
  }
  output.matrix <- cbind(output.matrix[, 1], temp)
  
  # otherwise we get duplicate model names in latex and html
  if (dots$output.type == 'ascii') {
    output.matrix <- rbind(c("", mod.names), output.matrix)
  }
  
  # add custom columns
  output.matrix <- customcolumns(output.matrix, custom.columns, custom.col.pos, 
      single.row = single.row, numcoef = nrow(m), groups = groups, 
      modelnames = TRUE)

  # attributes required for printing functions
  if ('include.attributes' %in% names(dots)) {
    if (dots$include.attributes) {
      attr(output.matrix, 'ci') <- ci
      attr(output.matrix, 'ci.test') <- ci.test
      attr(output.matrix, 'gof.names') <- gof.names
      attr(output.matrix, 'coef.names') <- coef.names
      attr(output.matrix, 'mod.names') <- mod.names
    }
  }
  
  return(output.matrix)
} 

# screenreg function
screenreg <- function(l,
                      file = NULL,
                      single.row = FALSE,
                      stars = c(0.001, 0.01, 0.05),
                      custom.model.names = NULL,
                      custom.coef.names = NULL,
                      custom.coef.map = NULL,
                      custom.gof.names = NULL,
                      custom.note = NULL,
                      digits = 2,
                      leading.zero = TRUE,
                      star.symbol = '*',
                      symbol = ".",
                      override.coef = 0,
                      override.se = 0,
                      override.pvalues = 0,
                      override.ci.low = 0, 
                      override.ci.up = 0,
                      omit.coef = NULL,
                      reorder.coef = NULL,
                      reorder.gof = NULL,
                      ci.force = "null",
                      ci.force.level = 0.95,
                      ci.test = 0,
                      groups = NULL,
                      custom.columns = NULL,
                      custom.col.pos = NULL,
                      column.spacing = 2,
                      outer.rule = "=",
                      inner.rule = "-",
                      ...) {
    
    # matrixreg produces the output matrix
    output.matrix <- matrixreg(l,
                               single.row = single.row,
                               stars = stars,
                               custom.model.names = custom.model.names,
                               custom.coef.names = custom.coef.names,
                               custom.coef.map = custom.coef.map,
                               custom.gof.names = custom.gof.names,
                               digits = digits,
                               leading.zero = leading.zero,
                               star.symbol = star.symbol,
                               symbol = symbol,
                               override.coef = override.coef,
                               override.se = override.se,
                               override.pvalues = override.pvalues,
                               override.ci.low = override.ci.low,
                               override.ci.up = override.ci.up,
                               omit.coef = omit.coef,
                               reorder.coef = reorder.coef,
                               reorder.gof = reorder.gof,
                               ci.force = ci.force,
                               ci.force.level = ci.force.level,
                               ci.test = ci.test,
                               groups = groups,
                               custom.columns = custom.columns,
                               custom.col.pos = custom.col.pos,
                               include.attributes = TRUE,
                               ...)
    
    gof.names <- attr(output.matrix, 'gof.names')
    coef.names <- attr(output.matrix, 'coef.names')
    mod.names <- attr(output.matrix, 'mod.names')
    ci <- attr(output.matrix, 'ci')
    ci.test <- attr(output.matrix, 'ci.test')
    
    # add spaces
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
    
    if (length(gof.names) > 0) {
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
    }
    
    # write table footer
    if (outer.rule != "") {
        string <- paste0(string, o.rule, "\n")
    }
    
    # stars note
    snote <- get_stars(pval = NULL,
                       stars = stars,
                       star.symbol = star.symbol,
                       symbol = symbol,
                       ci = ci,
                       ci.test = ci.test,
                       output = 'ascii')$note
    
    # custom note
    if (is.null(custom.note)) {
        note <- paste0(snote, "\n")
    } else if (custom.note == "") {
        note <- ""
    } else {
        note <- paste0(custom.note, "\n")
        note <- gsub("%stars", snote, note)
    }
    string <- paste0(string, note)
    
    #write to file
    if (is.null(file) || is.na(file)) {
        class(string) <- c("character", "texregTable")
        return(string)
    } else if (!is.character(file)) {
        stop("The 'file' argument must be a character string.")
    } else {
        sink(file)
        cat(string)
        sink()
        message(paste0("The table was written to the file '", file, "'.\n"))
    }
}


# texreg function
texreg <- function(l,
                   file = NULL,
                   single.row = FALSE,
                   stars = c(0.001, 0.01, 0.05),
                   custom.model.names = NULL,
                   custom.coef.names = NULL,
                   custom.coef.map = NULL,
                   custom.gof.names = NULL,
                   custom.note = NULL,
                   digits = 2,
                   leading.zero = TRUE,
                   symbol = "\\cdot",
                   override.coef = 0,
                   override.se = 0,
                   override.pvalues = 0,
                   override.ci.low = 0,
                   override.ci.up = 0,
                   omit.coef = NULL,
                   reorder.coef = NULL,
                   reorder.gof = NULL,
                   ci.force = "null",
                   ci.force.level = 0.95,
                   ci.test = 0,
                   groups = NULL,
                   custom.columns = NULL,
                   custom.col.pos = NULL,
                   bold = 0.00,
                   center = TRUE,
                   caption = "Statistical models",
                   caption.above = FALSE,
                   label = "table:coefficients",
                   booktabs = FALSE,
                   dcolumn = FALSE,
                   lyx = FALSE,
                   sideways = FALSE,
                   longtable = FALSE,
                   use.packages = TRUE,
                   table = TRUE,
                   no.margin = FALSE,
                   fontsize = NULL,
                   scalebox = NULL,
                   float.pos = "",
                   ...) {
  
  #check dcolumn vs. bold
  if (dcolumn == TRUE && bold > 0) {
    dcolumn <- FALSE
    msg <- paste("The dcolumn package and the bold argument cannot be used at", 
        "the same time. Switching off dcolumn.")
    if (length(stars) > 1 || stars == TRUE) {
      warning(paste(msg, "You should also consider setting stars = 0."))
    } else {
      warning(msg)
    }
  }
  
  # check longtable vs. sideways
  if (longtable == TRUE && sideways == TRUE) {
    sideways <- FALSE
    msg <- paste("The longtable package and sideways environment cannot be", 
        "used at the same time. You may want to use the pdflscape package.", 
        "Switching off sideways.")
    warning(msg)
  }
  
  # check longtable vs. float.pos
  if (longtable == TRUE && !(float.pos %in% c("", "l", "c", "r"))) {
    float.pos <- ""
    msg <- paste("When the longtable environment is used, the float.pos", 
        "argument can only take one of the \"l\", \"c\", \"r\", or \"\"", 
        "(empty) values. Setting float.pos = \"\".")
    warning(msg)
  }
  
  # check longtable vs. scalebox
  if (longtable == TRUE && !is.null(scalebox)) {
    scalebox <- NULL
    warning(paste("longtable and scalebox are not compatible. Setting", 
    "scalebox = NULL."))
  }

  # matrixreg produces the output matrix
  output.matrix <- matrixreg(l,
                             single.row = single.row,
                             stars = stars,
                             custom.model.names = custom.model.names,
                             custom.coef.names = custom.coef.names,
                             custom.coef.map = custom.coef.map,
                             custom.gof.names = custom.gof.names,
                             digits = digits,
                             leading.zero = leading.zero,
                             star.symbol = '*',
                             symbol = symbol,
                             override.coef = override.coef,
                             override.se = override.se,
                             override.pvalues = override.pvalues,
                             override.ci.low = override.ci.low, 
                             override.ci.up = override.ci.up,
                             omit.coef = omit.coef,
                             reorder.coef = reorder.coef,
                             reorder.gof = reorder.gof,
                             ci.force = ci.force,
                             ci.force.level = ci.force.level,
                             ci.test = ci.test,
                             groups = groups,
                             custom.columns = custom.columns,
                             custom.col.pos = custom.col.pos,
                             dcolumn = dcolumn,
                             bold = bold,
                             include.attributes = TRUE,
                             output.type = 'latex',
                             ...)
  
  gof.names <- attr(output.matrix, 'gof.names')
  coef.names <- attr(output.matrix, 'coef.names')
  mod.names <- attr(output.matrix, 'mod.names')
  ci <- attr(output.matrix, 'ci')
  ci.test <- attr(output.matrix, 'ci.test')
  
  # what is the optimal length of the labels?
  lab.list <- c(coef.names, gof.names)
  lab.length <- 0
  for (i in 1:length(lab.list)) {
    if (nchar(lab.list[i]) > lab.length) {
      lab.length <- nchar(lab.list[i])
    }
  }
  
  coltypes <- customcolumnnames(mod.names, custom.columns, custom.col.pos, 
      types = TRUE)
  mod.names <- customcolumnnames(mod.names, custom.columns, custom.col.pos, 
      types = FALSE)
  
  # define columns of the table (define now, add later)
  coldef <- ""
  if (no.margin == FALSE) {
    margin.arg <- ""
  } else {
    margin.arg <- "@{}"
  }
  coefcount <- 0
  for (i in 1:length(mod.names)) {
    if (coltypes[i] == "coef") {
      coefcount <- coefcount + 1
    }
    if (single.row == TRUE && coltypes[i] == "coef") {
      if (ci[coefcount] == FALSE) {
        separator <- ")"
      } else {
        separator <- "]"
      }
    } else {
      separator <- "."
    }
        if (coltypes[i] %in% c("coef", "customcol")) {
            alignmentletter <- "c"
        } else if (coltypes[i] == "coefnames") {
            alignmentletter <- "l"
        }
        if (dcolumn == FALSE) {
            coldef <- paste0(coldef, alignmentletter, margin.arg, " ")
        } else {
            if (coltypes[i] != "coef") {
                coldef <- paste0(coldef, alignmentletter, margin.arg, " ")
            } else {
                if (single.row == TRUE) {
                    dl <- compute.width(output.matrix[, i], left = TRUE, 
                                        single.row = TRUE, bracket = separator)
                    dr <- compute.width(output.matrix[, i], left = FALSE, 
                                        single.row = TRUE, bracket = separator)
                } else {
                    dl <- compute.width(output.matrix[, i], left = TRUE, 
                                        single.row = FALSE, bracket = separator)
                    dr <- compute.width(output.matrix[, i], left = FALSE, 
                                        single.row = FALSE, bracket = separator)
                }
                coldef <- paste0(coldef, "D{", separator, "}{", separator, "}{", 
                                 dl, separator, dr, "}", margin.arg, " ")
            }
        }
    }
    
    string <- "\n"
    linesep <- if (lyx) "\n\n" else "\n"
    
    # write table header
    if (use.packages == TRUE) {
        if (sideways == TRUE & table == TRUE) {
            string <- paste0(string, "\\usepackage{rotating}", linesep)
        }
        if (booktabs == TRUE) {
            string <- paste0(string, "\\usepackage{booktabs}", linesep)
        }
        if (dcolumn == TRUE) {
            string <- paste0(string, "\\usepackage{dcolumn}", linesep)
        }
        if (longtable == TRUE) {
            string <- paste0(string, "\\usepackage{longtable}", linesep)
        }
        if (dcolumn == TRUE || booktabs == TRUE || sideways == TRUE || 
            longtable == TRUE) {
            string <- paste0(string, linesep)
        }
    }
    
    if (longtable == TRUE) {
        if (center == TRUE) {
            string <- paste0(string, "\\begin{center}\n")
        }
        if (!is.null(fontsize)) {
            string <- paste0(string, "\\begin{", fontsize, "}", linesep)
        }
        if (float.pos == "") {
            string <- paste0(string, "\\begin{longtable}{", coldef, "}", linesep)
        } else {
            string <- paste0(string, "\\begin{longtable}[", float.pos, "]", linesep)
        }
    } else {  # table or sidewaystable
        if (table == TRUE) {
            if (sideways == TRUE) {
                t <- "sideways"
            } else {
                t <- ""
            }
            if (float.pos == "") {
                string <- paste0(string, "\\begin{", t, "table}", linesep)
            } else {
                string <- paste0(string, "\\begin{", t, "table}[", float.pos, "]", 
                                 linesep)
            }
            if (caption.above == TRUE) {
                string <- paste0(string, "\\caption{", caption, "}", linesep)
            }
            if (center == TRUE) {
                string <- paste0(string, "\\begin{center}", linesep)
            }
            if (!is.null(fontsize)) {
                string <- paste0(string, "\\begin{", fontsize, "}", linesep)
            }
            if (!is.null(scalebox)) {
                string <- paste0(string, "\\scalebox{", scalebox, "}{\n")
            }
        }
        string <- paste0(string, "\\begin{tabular}{", coldef, "}", linesep)
    }
    
    # horizontal rule above the table
    tablehead <- ""
    if (booktabs == TRUE) {
        tablehead <- paste0(tablehead, "\\toprule", linesep)
    } else {
        tablehead <- paste0(tablehead, "\\hline", linesep)
    }
    
    # specify model names
    tablehead <- paste0(tablehead, mod.names[1])
    if (dcolumn == TRUE) {
        for (i in 2:length(mod.names)) {
            if (coltypes[i] != "coef") {
                tablehead <- paste0(tablehead, " & ", mod.names[i])
            } else {
                tablehead <- paste0(tablehead, " & \\multicolumn{1}{c}{", mod.names[i], 
                                    "}")
            }
        }
    } else {
        for (i in 2:length(mod.names)) {
            tablehead <- paste0(tablehead, " & ", mod.names[i])
        }
    }
    
    # horizontal rule between model names and coefficients (define now, add later)
    if (booktabs == TRUE) {
        tablehead <- paste0(tablehead, " \\\\", linesep, "\\midrule", linesep)
    } else {
        tablehead <- paste0(tablehead, " \\\\", linesep, "\\hline", linesep)
    }
    if (longtable == FALSE) {
        string <- paste0(string, tablehead)
    }
    
    # stars note (define now, add later)
    snote <- get_stars(pval = NULL,
                       stars = stars,
                       star.symbol = '*',
                       symbol = symbol,
                       ci = ci,
                       ci.test = ci.test,
                       output = 'latex')$note
    
    if (is.null(fontsize)) {
        notesize <- "scriptsize"
    } else if (fontsize == "tiny" || fontsize == "scriptsize" || 
               fontsize == "footnotesize" || fontsize == "small") {
        notesize <- "tiny"
    } else if (fontsize == "normalsize") {
        notesize <- "scriptsize"
    } else if (fontsize == "large") {
        notesize <- "footnotesize"
    } else if (fontsize == "Large") {
        notesize <- "small"
    } else if (fontsize == "LARGE") {
        notesize <- "normalsize"
    } else if (fontsize == "huge") {
        notesize <- "large"
    } else if (fontsize == "Huge") {
        notesize <- "Large"
    }
    if (is.null(custom.note)) {
        if (snote == "") {
            note <- ""
        } else {
            note <- paste0("\\multicolumn{", length(mod.names), 
                           "}{l}{\\", notesize, "{", snote, "}}")
        }
    } else if (custom.note == "") {
        note <- ""
    } else {
        note <- paste0("\\multicolumn{", length(mod.names), 
                       "}{l}{\\", notesize, "{", custom.note, "}}")
        note <- gsub("%stars", snote, note, perl = TRUE)
    }
    if (note != "") {
        if (longtable == TRUE) {  # longtable requires line break after note/caption
            note <- paste0(note, "\\\\", linesep)
        } else {
            note <- paste0(note, linesep)
        }
    }
    
    # bottom rule (define now, add later)
    if (booktabs == TRUE) {
        bottomline <- paste0("\\bottomrule", linesep)
    } else {
        bottomline <- paste0("\\hline", linesep)
    }
    
    # write table header (and footer, in the case of longtable)
    if (longtable == TRUE) {
        if (caption.above == TRUE) {
            string <- paste0(string, "\\caption{", caption, "}", linesep, "\\label{", 
                             label, "}\\\\", linesep, tablehead, "\\endfirsthead", linesep, 
                             tablehead, "\\endhead", linesep, bottomline, "\\endfoot", linesep, 
                             bottomline, note, "\\endlastfoot", linesep)
        } else {
            string <- paste0(string, tablehead, "\\endfirsthead", linesep, tablehead, 
                             "\\endhead", linesep, bottomline, "\\endfoot", linesep, bottomline, 
                             note, "\\caption{", caption, "}", linesep, "\\label{", label, "}", 
                             linesep, "\\endlastfoot \\\\", linesep)
        }
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
                string <- paste0(string, " \\\\", linesep)
            } else {
                string <- paste0(string, " & ")
            }
        }
    }
    
    if (length(gof.names) > 0) {
        # lower mid rule
        if (booktabs == TRUE) {
            string <- paste0(string, "\\midrule", linesep)
        } else {
            string <- paste0(string, "\\hline", linesep)
        }
        
        # write GOF block
        for (i in (length(output.matrix[, 1]) - (length(gof.names) - 1)):
             (length(output.matrix[, 1]))) {
            for (j in 1:length(output.matrix[1, ])) {
                string <- paste0(string, output.matrix[i, j])
                if (j == length(output.matrix[1, ])) {
                    string <- paste0(string, " \\\\", linesep)
                } else {
                    string <- paste0(string, " & ")
                }
            }
        }
    }
    
    # write table footer
    if (longtable == FALSE) {
        string <- paste0(string, bottomline)
        string <- paste0(string, note, "\\end{tabular}", linesep)
    }
    
    # take care of center, scalebox and table environment
    if (longtable == TRUE) {
        string <- paste0(string, "\\end{longtable}", linesep)
        if (!is.null(fontsize)) {
            string <- paste0(string, "\\end{", fontsize, "}", linesep)
        }
        if (center == TRUE) {
            string <- paste0(string, "\\end{center}", linesep)
        }
    } else if (table == TRUE) {
        if (!is.null(fontsize)) {
            string <- paste0(string, "\\end{", fontsize, "}", linesep)
        }
        if (!is.null(scalebox)) {
            string <- paste0(string, "}", linesep)
        }
        if (caption.above == FALSE) {
            string <- paste0(string, "\\caption{", caption, "}", linesep)
        }
        string <- paste0(string, "\\label{", label, "}", linesep)
        if (center == TRUE) {
            string <- paste0(string, "\\end{center}", linesep)
        }
        if (sideways == TRUE) {
            t <- "sideways"
        } else {
            t <- ""
        }
        string <- paste0(string, "\\end{", t, "table}", linesep)
    }
    
    if (is.null(file) || is.na(file)) {
        class(string) <- c("character", "texregTable")
        return(string)
    } else if (!is.character(file)) {
        stop("The 'file' argument must be a character string.")
    } else {
        sink(file)
        cat(string)
        sink()
        message(paste0("The table was written to the file '", file, "'.\n"))
    }
}


# htmlreg function
htmlreg <- function(l,
                    file = NULL,
                    single.row = FALSE,
                    stars = c(0.001, 0.01, 0.05),
                    custom.model.names = NULL, 
                    custom.coef.names = NULL,
                    custom.coef.map = NULL,
                    custom.gof.names = NULL,
                    custom.note = NULL,
                    digits = 2,
                    leading.zero = TRUE,
                    star.symbol = '*',
                    symbol = "&middot;",
                    override.coef = 0,
                    override.se = 0,
                    override.pvalues = 0,
                    override.ci.low = 0,
                    override.ci.up = 0,
                    omit.coef = NULL,
                    reorder.coef = NULL,
                    reorder.gof = NULL,
                    ci.force = "null",
                    ci.force.level = 0.95,
                    ci.test = 0,
                    groups = NULL,
                    custom.columns = NULL,
                    custom.col.pos = NULL,
                    bold = 0.00,
                    center = TRUE,
                    caption = "Statistical models",
                    caption.above = FALSE,
                    inline.css = TRUE,
                    doctype = TRUE,
                    html.tag = FALSE,
                    head.tag = FALSE,
                    body.tag = FALSE,
                    indentation = "",
                    vertical.align.px = 0,
                    ...) {
    
    # inline CSS definitions
    if (inline.css == TRUE) {
        css.table <- " style=\"border: none;\""
        css.th <- paste0(" style=\"text-align: left; border-top: 2px solid ", 
                         "black; border-bottom: 1px solid black; padding-right: 12px;\"")
        css.midrule <- " style=\"border-top: 1px solid black;\""
        css.bottomrule <- " style=\"border-bottom: 2px solid black;\""
        css.bottomrule.nogof <- paste(" style=\"padding-right: 12px;",
                                      "border-bottom: 2px solid black;\"")
        css.td <- " style=\"padding-right: 12px; border: none;\""
        css.caption <- ""
        css.sup <- paste0(" style=\"vertical-align: ", vertical.align.px, "px;\"")
    } else {
        css.table <- ""
        css.th <- ""
        css.midrule <- ""
        css.bottomrule <- ""
        css.td <- ""
        css.caption <- ""
        css.sup <- ""
    }
    
    # matrixreg produces the output matrix
    output.matrix <- matrixreg(l,
                               single.row = single.row, 
                               stars = stars,
                               custom.model.names = custom.model.names,
                               custom.coef.names = custom.coef.names,
                               custom.coef.map = custom.coef.map,
                               custom.gof.names = custom.gof.names,
                               digits = digits,
                               leading.zero = leading.zero,
                               star.symbol = star.symbol,
                               symbol = symbol,
                               override.coef = override.coef,
                               override.se = override.se,
                               override.pvalues = override.pvalues,
                               override.ci.low = override.ci.low, 
                               override.ci.up = override.ci.up,
                               omit.coef = omit.coef,
                               reorder.coef = reorder.coef,
                               reorder.gof = reorder.gof,
                               ci.force = ci.force,
                               ci.force.level = ci.force.level,
                               ci.test = ci.test,
                               groups = groups,
                               custom.columns = custom.columns,
                               custom.col.pos = custom.col.pos,
                               bold = bold,
                               include.attributes = TRUE,
                               output.type = 'html',
                               css.sup = css.sup,
                               ...)
    
    gof.names <- attr(output.matrix, 'gof.names')
    coef.names <- attr(output.matrix, 'coef.names')
    mod.names <- attr(output.matrix, 'mod.names')
    ci <- attr(output.matrix, 'ci')
    ci.test <- attr(output.matrix, 'ci.test')
    
    coltypes <- customcolumnnames(mod.names, custom.columns, custom.col.pos, 
                                  types = TRUE)
    mod.names <- customcolumnnames(mod.names, custom.columns, custom.col.pos, 
                                   types = FALSE)
    
    
    # write table header
    if (single.row == TRUE) {
        numcols <- 2 * length(mod.names)
    } else {
        numcols <- length(mod.names)
    }
    
    if (doctype == TRUE) {
        doct <- paste0("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 ", 
                       "Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n")
    } else {
        doct <- ""
    }
    
    # determine indentation for table
    ind <- indentation
    if (html.tag == TRUE) {
        h.ind <- indentation
    } else {
        h.ind <- ""
    }
    if (body.tag == TRUE) {
        b.ind <- indentation
    } else {
        b.ind <- ""
    }
    if (head.tag == TRUE) {
        d.ind <- indentation
    } else {
        d.ind <- ""
    }
    
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
            h.ind, d.ind, ind, ind, "vertical-align: ", vertical.align.px, "px;\n", 
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
        h.ind, b.ind, ind, "<tr>\n"
    )
    
    # specify model names (header row)
    for (i in 1:length(mod.names)) {
        string <- paste0(string, 
                         h.ind, b.ind, ind, ind, "<th", css.th, "><b>", mod.names[i], 
                         "</b></th>\n")
    }
    string <- paste0(string, h.ind, b.ind, ind, "</tr>\n")
    
    # write coefficients to string object
    coef.length <- length(output.matrix[, 1]) - length(gof.names)
    for (i in 1:coef.length) {
        string <- paste0(string, h.ind, b.ind, ind, "<tr>\n")
        for (j in 1:length(output.matrix[1, ])) {
            if (length(gof.names) == 0 && i == coef.length) { # no GOF block
                if (inline.css == TRUE) {
                    br <- css.bottomrule.nogof
                } else {
                    br <- " class=\"bottomRule\""
                }
                string <- paste0(string, h.ind, b.ind, ind, ind, "<td", br, ">", 
                                 output.matrix[i,j], "</td>\n")
            } else { # GOF block present
                string <- paste0(string, h.ind, b.ind, ind, ind, "<td", css.td, ">", 
                                 output.matrix[i,j], "</td>\n")
            }
        }
        string <- paste0(string, h.ind, b.ind, ind, "</tr>\n")
    }
    
    if (length(gof.names) > 0) {
        # write GOF block
        for (i in (length(output.matrix[, 1]) - (length(gof.names) - 1)):
             (length(output.matrix[, 1]))) {
            string <- paste0(string, h.ind, b.ind, ind, "<tr>\n")
            for (j in 1:length(output.matrix[1, ])) {
                if (i == length(output.matrix[, 1]) - (length(gof.names) - 1)) {
                    if (inline.css == TRUE) {
                        mr <- css.midrule
                    } else {
                        mr <- " class=\"midRule\""  # add mid rule via style sheets
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
    }
    
    # stars note
    snote <- get_stars(pval = NULL,
                       stars = stars,
                       star.symbol = star.symbol,
                       symbol = symbol,
                       ci = ci,
                       ci.test = ci.test,
                       css.sup = css.sup,
                       output = 'html')$note
    
    if (is.null(custom.note)) {
        note <- snote
    } else if (custom.note == "") {
        note <- ""
    } else {
        note <- custom.note
        note <- gsub("%stars", snote, note)
    }
    if (note != "") {
        string <- paste0(string,
                         h.ind,
                         b.ind,
                         ind,
                         "<tr>\n",
                         h.ind,
                         b.ind,
                         ind,
                         ind,
                         "<td",
                         css.td,
                         " colspan=\"",
                         (1 + length(mod.names)),
                         "\"><span style=\"font-size:0.8em\">",
                         note,
                         "</span></td>\n",
                         h.ind,
                         b.ind,
                         ind,
                         "</tr>\n")
    }
    
    # write table footer
    string <- paste0(string, h.ind, b.ind, "</table>\n")
    if (body.tag == TRUE) {
        string <- paste0(string, h.ind, "</body>\n")
    }
    if (html.tag == TRUE) {
        string <- paste0(string, "</html>\n")
    }
    
    if (is.null(file) || is.na(file)) {
        class(string) <- c("character", "texregTable")
        return(string)
    } else if (!is.character(file)) {
        stop("The 'file' argument must be a character string.")
    } else {
        sink(file)
        cat(string)
        sink()
        message(paste0("The table was written to the file '", file, "'.\n"))
    }
}

# csvreg function
csvreg <- function(l,
                   file,
                   stars = c(0.001, 0.01, 0.05),
                   custom.model.names = NULL,
                   custom.coef.names = NULL,
                   custom.coef.map = NULL,
                   custom.gof.names = NULL,
                   custom.note = NULL,
                   digits = 2,
                   leading.zero = TRUE,
                   star.symbol = '*',
                   symbol = ".",
                   override.coef = 0,
                   override.se = 0,
                   override.pvalues = 0,
                   override.ci.low = 0,
                   override.ci.up = 0,
                   omit.coef = NULL,
                   reorder.coef = NULL,
                   reorder.gof = NULL,
                   ci.force = "null",
                   ci.force.level = 0.95,
                   ci.test = 0,
                   groups = NULL,
                   custom.columns = NULL,
                   custom.col.pos = NULL,
                   caption = 'Statistical Models',
                   ...) {
    
    # matrixreg produces the output matrix
    output.matrix <- matrixreg(l,
                               stars = stars,
                               custom.model.names = custom.model.names,
                               custom.coef.names = custom.coef.names,
                               custom.coef.map = custom.coef.map,
                               custom.gof.names = custom.gof.names,
                               digits = digits,
                               leading.zero = leading.zero,
                               star.symbol = star.symbol,
                               symbol = symbol, 
                               override.coef = override.coef,
                               override.se = override.se,
                               override.pvalues = override.pvalues,
                               override.ci.low = override.ci.low,
                               override.ci.up = override.ci.up,
                               omit.coef = omit.coef,
                               reorder.coef = reorder.coef,
                               reorder.gof = reorder.gof,
                               ci.force = ci.force,
                               ci.force.level = ci.force.level,
                               ci.test = ci.test,
                               groups = groups,
                               custom.columns = custom.columns,
                               custom.col.pos = custom.col.pos,
                               include.attributes = TRUE,
                               ...)
    
    # attributes
    ci <- attr(output.matrix, 'ci')
    ci.test <- attr(output.matrix, 'ci.test')
    
    # append notes to bottom of table 
    out <- output.matrix
    if (is.character(caption) && (caption != '')) {
        out <- rbind(out, c('Caption: ', caption, rep('', ncol(output.matrix) - 2)))
    }
    snote <- get_stars(pval = NULL,
                       stars = stars,
                       star.symbol = star.symbol,
                       symbol = symbol,
                       ci = ci,
                       ci.test = ci.test,
                       output = 'ascii')$note
    if (trimws(snote) != '') {
        out <- rbind(out, c('Note: ', snote, rep('', ncol(output.matrix) - 2)))
    } 
    if (is.character(custom.note) && (custom.note != '')) {
        out <- rbind(out, c('Note: ', custom.note, rep('', ncol(output.matrix) - 2)))
    }
    out <- as.data.frame(out)
    
    # write csv to file
    if (!is.character(file)) {
        stop('file must be a character string')
    } else {
        write.table(out,
                    file = file,
                    sep = ',',
                    quote = TRUE,
                    col.names = FALSE,
                    row.names = FALSE)
    }
}

# Microsoft Word function
wordreg <- function(l,
                    file = NULL,
                    single.row = FALSE,
                    stars = c(0.001, 0.01, 0.05),
                    custom.model.names = NULL,
                    custom.coef.names = NULL,
                    custom.coef.map = NULL,
                    custom.gof.names = NULL,
                    digits = 2,
                    leading.zero = TRUE,
                    star.symbol = star.symbol,
                    symbol = ".",
                    override.coef = 0,
                    override.se = 0,
                    override.pvalues = 0,
                    override.ci.low = 0,
                    override.ci.up = 0,
                    omit.coef = NULL,
                    reorder.coef = NULL,
                    reorder.gof = NULL,
                    ci.force = "null",
                    ci.force.level = 0.95,
                    ci.test = 0,
                    groups = NULL,
                    custom.columns = NULL,
                    custom.col.pos = NULL,
                    ...) {
    
    if (!'rmarkdown' %in% row.names(installed.packages())) {
        stop(paste("The wordreg function requires the 'rmarkdown' package.",
                   "Install it and try again."))
    }
    if (is.null(file)) {
        stop("'file' must be a valid file path.")
    }
    mat <- matrixreg(l, 
                     single.row = single.row,
                     stars = stars,
                     custom.model.names = custom.model.names,
                     custom.coef.names = custom.coef.names, 
                     custom.coef.map = custom.coef.map,
                     custom.gof.names = custom.gof.names,
                     digits = digits,
                     leading.zero = leading.zero,
                     star.symbol = '*', # produces error if fed star.symbol itself
                     symbol = symbol,
                     override.coef = override.coef,
                     override.se = override.se,
                     override.pvalues = override.pvalues,
                     override.ci.low = override.ci.low,
                     override.ci.up = override.ci.up,
                     omit.coef = omit.coef,
                     reorder.coef = reorder.coef,
                     reorder.gof = reorder.gof,
                     ci.force = ci.force,
                     ci.force.level = ci.force.level,
                     ci.test = ci.test,
                     groups = groups,
                     custom.columns = custom.columns,
                     custom.col.pos = custom.col.pos,
                     output.type = 'ascii',
                     include.attributes = FALSE
    )
    wd <- getwd()
    f = tempfile(fileext = '.Rmd')
    cat(file = f, '```{r, echo = FALSE}
                    knitr::kable(mat)
                    ```', append = TRUE)
    rmarkdown::render(f, output_file = paste0(wd, "/", file))
}
