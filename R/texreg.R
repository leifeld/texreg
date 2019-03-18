# The texreg package was written by Philip Leifeld.
# Please use the issue tracker at http://github.com/leifeld/texreg
# for bug reports, help, or feature requests.


# texreg functions -------------------------------------------------------------

# csvreg function
csvreg <- function(l,
                   file,
                   stars = c(0.001, 0.01, 0.05),
                   custom.model.names = NULL,
                   custom.coef.names = NULL,
                   custom.coef.map = NULL,
                   custom.gof.names = NULL,
                   custom.gof.rows = NULL,
                   custom.note = NULL,
                   digits = 2,
                   leading.zero = TRUE,
                   star.symbol = "*",
                   symbol = ".",
                   override.coef = 0,
                   override.se = 0,
                   override.pvalues = 0,
                   override.ci.low = 0,
                   override.ci.up = 0,
                   omit.coef = NULL,
                   reorder.coef = NULL,
                   reorder.gof = NULL,
                   ci.force = FALSE,
                   ci.force.level = 0.95,
                   ci.test = 0,
                   groups = NULL,
                   custom.columns = NULL,
                   custom.col.pos = NULL,
                   caption = "Statistical Models",
                   ...) {

  # matrixreg produces the output matrix
  output.matrix <- matrixreg(l,
                             stars = stars,
                             custom.model.names = custom.model.names,
                             custom.coef.names = custom.coef.names,
                             custom.coef.map = custom.coef.map,
                             custom.gof.names = custom.gof.names,
                             custom.gof.rows = custom.gof.rows,
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
                             trim = TRUE,
                             ...)

  # attributes
  ci <- attr(output.matrix, "ci")
  ci.test <- attr(output.matrix, "ci.test")

  # append notes to bottom of table
  out <- output.matrix
  if (is.character(caption) && (caption != "")) {
    out <- rbind(out, c("Caption: ", caption, rep("", ncol(output.matrix) - 2)))
  }
  snote <- get_stars(pval = NULL,
                     stars = stars,
                     star.symbol = star.symbol,
                     symbol = symbol,
                     ci = ci,
                     ci.test = ci.test,
                     output = "ascii")$note
  if (trimws(snote) != "") {
    out <- rbind(out, c("Note: ", snote, rep("", ncol(output.matrix) - 2)))
  }
  if (is.character(custom.note) && (custom.note != "")) {
    out <- rbind(out, c("Note: ", custom.note, rep("", ncol(output.matrix) - 2)))
  }
  out <- as.data.frame(out)

  # write csv to file
  if (!is.character(file)) {
    stop("File must be a character string.")
  } else {
    write.table(out,
                file = file,
                sep = ",",
                quote = TRUE,
                col.names = FALSE,
                row.names = FALSE)
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
                    custom.gof.rows = NULL,
                    custom.note = NULL,
                    digits = 2,
                    leading.zero = TRUE,
                    star.symbol = "*",
                    symbol = "&middot;",
                    override.coef = 0,
                    override.se = 0,
                    override.pvalues = 0,
                    override.ci.low = 0,
                    override.ci.up = 0,
                    omit.coef = NULL,
                    reorder.coef = NULL,
                    reorder.gof = NULL,
                    ci.force = FALSE,
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
                             custom.gof.rows = custom.gof.rows,
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
                             trim = FALSE,
                             output.type = "html",
                             css.sup = css.sup,
                             ...)

  gof.names <- attr(output.matrix, "gof.names")
  coef.names <- attr(output.matrix, "coef.names")
  mod.names <- attr(output.matrix, "mod.names")
  ci <- attr(output.matrix, "ci")
  ci.test <- attr(output.matrix, "ci.test")

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
                     output = "html")$note

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

# huxtablereg function
huxtablereg <- function(l,
                        single.row = FALSE,
                        stars = c(0.001, 0.01, 0.05),
                        custom.model.names = NULL,
                        custom.coef.names = NULL,
                        custom.coef.map = NULL,
                        custom.gof.names = NULL,
                        custom.gof.rows = NULL,
                        digits = 2,
                        leading.zero = TRUE,
                        star.symbol = star.symbol,
                        symbol = "+",
                        override.coef = 0,
                        override.se = 0,
                        override.pvalues = 0,
                        override.ci.low = 0,
                        override.ci.up = 0,
                        omit.coef = NULL,
                        reorder.coef = NULL,
                        reorder.gof = NULL,
                        ci.force = FALSE,
                        ci.force.level = 0.95,
                        ci.test = 0,
                        groups = NULL,
                        custom.columns = NULL,
                        custom.col.pos = NULL,
                        ...)  {
  if (!requireNamespace("huxtable", quietly = TRUE)) {
    stop("huxtablereg requires the 'huxtable' package to be installed.\n",
         "To do this, enter 'install.packages(\"huxtable\")'.")
  }

  mr.call <- match.call(expand.dots = FALSE)
  mr.call[[1L]] <- quote(texreg::matrixreg)
  mr.call$include.attributes <- TRUE
  mr.call$trim <- TRUE
  mx <- eval(mr.call)

  gof.names <- attr(mx, "gof.names")
  coef.names <- attr(mx, "coef.names")
  ci <- attr(mx, "ci")

  hx <- huxtable::as_hux(mx, add_colnames = FALSE, autoformat = TRUE)
  huxtable::align(hx)[-1, -1] <- "right"
  coef.rows <- which(as.matrix(hx[, 1]) %in% coef.names)
  hx <- huxtable::set_align(hx, coef.rows, -1, ".")
  if (!single.row) {
    hx <- huxtable::set_align(hx, coef.rows + 1, c(FALSE, !ci), ".")
  }
  hx <- huxtable::set_number_format(hx, coef.rows, -1, digits)

  gof.rows <- as.matrix(hx[, 1]) %in% gof.names
  hx <- huxtable::set_align(hx, gof.rows, -1, ".")

  return(hx)
}

# matrixreg function
matrixreg <- function(l,
                      single.row = FALSE,
                      stars = c(0.001, 0.01, 0.05),
                      custom.model.names = NULL,
                      custom.coef.names = NULL,
                      custom.coef.map = NULL,
                      custom.gof.names = NULL,
                      custom.gof.rows = NULL,
                      digits = 2,
                      leading.zero = TRUE,
                      star.symbol = "*",
                      symbol = ".",
                      override.coef = 0,
                      override.se = 0,
                      override.pvalues = 0,
                      override.ci.low = 0,
                      override.ci.up = 0,
                      omit.coef = NULL,
                      reorder.coef = NULL,
                      reorder.gof = NULL,
                      ci.force = FALSE,
                      ci.force.level = 0.95,
                      ci.test = 0,
                      bold = 0,
                      groups = NULL,
                      custom.columns = NULL,
                      custom.col.pos = NULL,
                      dcolumn = TRUE,
                      output.type = c("ascii", "latex", "html"),
                      include.attributes = FALSE,
                      trim = FALSE,
                      ...) {

  # unnamed arguments to environment
  dots <- list(...)

  # extract coefficients, SEs, GOFs etc. from the models and save in a list called 'models'
  models <- get.data(l, ...)

  # override coefs, SEs, p-values, and/or CIs if provided
  models <- override(models = models,
                     override.coef = override.coef,
                     override.se = override.se,
                     override.pvalues = override.pvalues,
                     override.ci.low = override.ci.low,
                     override.ci.up = override.ci.up)

  # convert LaTeX markup to text
  if (output.type[1] != "latex") {
    for (i in 1:length(models)) {
      # replace markup in GOF names
      if (output.type[1] == "html") {
        r <- paste0("<sup", dots$css.sup, ">2</sup>")
      } else if (output.type[1] == "ascii") {
        r <- "^2"
      } else {
        stop("'output.type' must be 'latex', 'html', or 'ascii'.")
      }
      models[[i]]@gof.names <- gsub("\\$\\^2\\$", r, models[[i]]@gof.names)
      models[[i]]@gof.names <- gsub("\\\\ ", " ", models[[i]]@gof.names)
      models[[i]]@gof.names <- gsub("\\ ", " ", models[[i]]@gof.names)

      # replace Greek letters in coefficient names
      models[[i]]@coef.names <- gsub("\\$\\\\rho\\$", "rho", models[[i]]@coef.names)
      models[[i]]@coef.names <- gsub("\\$\\\\lambda\\$", "lambda", models[[i]]@coef.names)
      models[[i]]@coef.names <- gsub("\\$\\\\mu\\$", "mu", models[[i]]@coef.names)
      models[[i]]@coef.names <- gsub("\\$\\\\nu\\$", "nu", models[[i]]@coef.names)
      models[[i]]@coef.names <- gsub("\\$\\\\tau\\$", "tau", models[[i]]@coef.names)
      models[[i]]@coef.names <- gsub("\\$\\\\sigma\\$", "sigma", models[[i]]@coef.names)
    }
  }

  # create confidence intervals using ci.force argument if necessary
  if (is.logical(ci.force) && length(ci.force) == 1) {
    ci.force <- rep(ci.force, length(models))
  }
  if (!is.logical(ci.force)) {
    stop("The 'ci.force' argument must be a vector of logical values.")
  }
  if (length(ci.force) != length(models)) {
    stop(paste("There are", length(models), "models and", length(ci.force), "ci.force values."))
  }
  if (is.null(ci.force.level) ||
      length(ci.force.level) != 1 ||
      !is.numeric(ci.force.level) ||
      ci.force.level > 1 ||
      ci.force.level < 0) {
    stop("'ci.force.level' must be a single value between 0 and 1.")
  }
  for (i in 1:length(models)) {
    if (ci.force[i] == TRUE && length(models[[i]]@se) > 0) {
      z <- qnorm(1 - ((1 - ci.force.level) / 2))
      upper <- models[[i]]@coef + (z * models[[i]]@se)
      lower <- models[[i]]@coef - (z * models[[i]]@se)
      models[[i]]@ci.low <- lower
      models[[i]]@ci.up <- upper
      models[[i]]@se <- numeric(0)
      models[[i]]@pvalues <- numeric(0)
    }
  }

  # extract names of the goodness-of-fit statistics (before adding any GOF rows)
  gof.names <- character()  # GOF names of all models in one vector
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

  # correct duplicate coefficient names (add " (1)", " (2)" etc.)
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
          models[[i]]@coef.names[k] <- paste0(models[[i]]@coef.names[k], " (", count, ")")
          count <- count + 1
        }
      }
    }
  }

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
      coef_i <- cbind(cf, cil, ciu)
    } else {
      if (length(se) > 0 && length(pv) > 0) {
        coef_i <- cbind(cf, se, pv)
      } else if (length(se) > 0 && length(pv) == 0) {
        # p-values not provided -> use p-values of 0.99
        coef_i <- cbind(cf, se, rep(0.99, length(cf)))
      } else if (length(se) == 0 && length(pv) > 0) {
        coef_i <- cbind(cf, rep(NA, length(cf)), pv)
      } else {
        # not even SEs provided
        coef_i <- cbind(cf, rep(NA, length(cf)), rep(0.99, length(cf)))
      }
    }
    rownames(coef_i) <- models[[i]]@coef.names
    coefs[[i]] <- coef_i
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

  # merge the coefficient tables into a new table called 'm'
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
    for (j in 1:ncol(m)) {
      m.temp[new.row, j] <- m[i, j]
    }
  }
  rownames(m.temp) <- coef.order
  colnames(m.temp) <- colnames(m)
  m <- m.temp
  rm(m.temp)

  # replace GOF names by custom names
  if (is.null(custom.gof.names)) {
    # do nothing
  } else if (!is.character(custom.gof.names)) {
    stop("Custom GOF names must be provided as a character vector.")
  } else if (length(custom.gof.names) != length(gof.names)) {
    stop(paste("There are", length(gof.names),
               "GOF statistics, but you provided", length(custom.gof.names),
               "custom names for them."))
  } else {
    custom.gof.names[is.na(custom.gof.names)] <- rownames(gofs)[is.na(custom.gof.names)]
    rownames(gofs) <- custom.gof.names
  }

  # add row names as first column to GOF block and format values as character strings
  if (dcolumn == FALSE && output.type[1] == "latex") {
    dollar <- "$"
  } else {
    dollar <- ""
  }
  gof.matrix <- matrix(nrow = nrow(gofs), ncol = ncol(gofs) + 1)  # including labels
  if (nrow(gof.matrix) > 0) {
    for (i in 1:nrow(gofs)) {
      # replace symbols in latex
      if (output.type[1] == "latex") {
        gof.matrix[i, 1] <- rownames(replaceSymbols(gofs))[i]
      } else {
        gof.matrix[i, 1] <- rownames(gofs)[i]
      }
      for (j in 1:ncol(gofs)) {
        strg <- coeftostring(gofs[i, j], leading.zero, digits = decimal.matrix[i, j])
        gof.matrix[i, j + 1] <- paste0(dollar, strg, dollar)
      }
    }
  }

  # add custom GOF rows
  if (!is.null(custom.gof.rows) && !is.na(custom.gof.rows)) {
    if (class(custom.gof.rows) != "list") {
      stop("The 'custom.gof.rows' argument is ignored because it is not a list.")
    }
    for (i in length(custom.gof.rows):1) {
      if (length(custom.gof.rows[[i]]) != ncol(gofs)) {
        stop("Custom GOF row ", i, " has a different number of values than there are models.")
      } else {
        if (!is.numeric(custom.gof.rows[[i]]) && !is.character(custom.gof.rows[[i]])) {
          custom.gof.rows[[i]] <- as.character(custom.gof.rows[[i]]) # cast into character if unknown class, such as factor or logical
        }
        # auto-detect decimal setting
        if (!is.numeric(custom.gof.rows)) { # put NA for character objects, for example fixed effects
          dec <- NA
        } else if (all(custom.gof.rows %% 1 == 0)) { # put 0 for integers
          dec <- 0
        } else { # put the respective decimal places if numeric but not integer
          dec <- digits
        }
        newValues <- sapply(custom.gof.rows[[i]], function(x) { # format the different values of the new row
          if (is.character(x) && output.type[1] == "latex" && isTRUE(dcolumn)) {
            paste0("\\multicolumn{1}{c}{", x, "}")
          } else {
            paste0(dollar, coeftostring(x, leading.zero, digits = dec), dollar)
          }
        })
        gof.matrix <- rbind(c(names(custom.gof.rows)[i], newValues), gof.matrix) # insert GOF name and GOF values as new row
      }
    }
  }

  # apply custom coefficient map using 'custom.coef.map' argument
  if (!is.null(custom.coef.map)) {
    # sanity checks
    if (class(custom.coef.map) != "list" || is.null(names(custom.coef.map))) {
      stop("'custom.coef.map' must be a named list.")
    }
    if (!any(names(custom.coef.map) %in% row.names(m))) {
      stop(paste("None of the coefficient names supplied in 'custom.coef.map'",
                 "appear to be in your models."))
    }

    # when user supplies NA as destination, replace with origin
    idx <- is.na(custom.coef.map)
    custom.coef.map[idx] <- names(custom.coef.map)[idx]

    # subset of coefficients to keep
    origin <- names(custom.coef.map)[names(custom.coef.map) %in% row.names(m)]
    destination <- unlist(custom.coef.map[origin])
    out <- m[origin, , drop = FALSE] # drop: otherwise R converts to numeric if a single coefficient is passed

    # rename
    row.names(out) <- destination
  } else { # use 'omit.coef' and 'custom.coef.names' if available
    # omit
    if (!is.null(omit.coef)) {
      if (!is.character(omit.coef) || is.na(omit.coef)) {
        stop("'omit.coef' must be a character object.")
      }
      idx <- !grepl(omit.coef, row.names(m), perl = TRUE)
      if (all(!idx)) {
        stop("You tried to remove all coefficients using 'omit.coef'.")
      }
    } else {
      idx <- rep(TRUE, nrow(m))
    }

    # rename
    if (!is.null(custom.coef.names)) {
      if (!is.character(custom.coef.names)) {
        stop("'custom.coef.names' must be a character vector.")
      }
      if (!length(custom.coef.names) %in% c(nrow(m), sum(idx))) { # check length
        if (nrow(m) == sum(idx)) {
          stop("'custom.coef.names' must be a character vector of length ", nrow(m), ".")
        } else {
          stop("'custom.coef.names' must be a character vector of length ", sum(idx), " or ", nrow(m), ".")
        }
      }
      # user submits number of custom names after omission
      if (length(custom.coef.names) == sum(idx)) {
        custom.coef.names <- custom.coef.names
      } else { # user submits number of custom names before omission
        custom.coef.names <- custom.coef.names[idx]
      }
    } else {
      custom.coef.names <- row.names(m)[idx]
    }

    # output
    m <- m[idx, , drop = FALSE]
    row.names(m) <- custom.coef.names
  }
  m <- as.data.frame(m)

  # reorder GOF block using reorder.gof argument
  gof.matrix <- reorder(gof.matrix, reorder.gof)

  # resort matrix and conflate duplicate entries:
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
  unique.names <- unique(rownames(m))               # unique row names in m
  num.unique <- length(unique.names)                # count these unique names
  orig.width <- ncol(m)                             # number of columns in m
  q <- matrix(nrow = 0, ncol = orig.width)          # new matrix with same width
  for (i in 1:num.unique) {                         # go through unique row names
    rows <- matrix(NA, nrow = 0, ncol = orig.width) # create matrix where re-
    # arranged rows will be stored
    for (j in 1:orig.width) {                       # go through columns in m
      current.name <- unique.names[i]               # save row name
      nonNa <- m[rownames(m) == current.name, j]    # create a vector of values
      # with same rowname in the col
      nonNa <- nonNa[!is.na(nonNa)]                 # retain only non-NA values
      for (k in 1:length(nonNa)) {                  # go through non-NA values
        if (k > dim(rows)[1]) {                     # add an NA-only row in which
          rows <- rbind(rows, rep(NA, orig.width))  # the values are stored
          rownames(rows)[k] <- unique.names[i]      # also add the row name
        }
        rows[k, j] <- nonNa[k]                      # actually store the value
      }
    }
    q <- rbind(q, rows)                             # add the new row(s) to q
  }
  m <- q

  # decide if default or custom model names should be used
  if (!is.null(custom.model.names) && !is.character(custom.model.names)) {
    stop("Model names in 'custom.model.names' must be specified as a character vector.")
  } else if (!is.null(custom.model.names) && length(custom.model.names) != length(models)) {
    stop(paste("There are", length(models), "models, but you provided",
               length(custom.model.names), "name(s) for them."))
  }
  mod.names <- character(length(models))
  for (i in 1:length(mod.names)) {
    if (!is.null(custom.model.names) && !is.na(custom.model.names[i]) && custom.model.names[i] != "") {
      mod.names[i] <- custom.model.names[i]
    } else if (!is.null(names(l)) && !is.na(names(l)[i]) && names(l)[i] != "") {
      mod.names[i] <- names(l)[i]
    } else if (!is.null(models[[i]]@model.name) &&
               !is.na(models[[i]]@model.name) &&
               !length(models[[i]]@model.name) == 0 &&
               models[[i]]@model.name != "") {
      mod.names[i] <- models[[i]]@model.name
    } else {
      mod.names[i] <- paste("Model", i)
    }
  }

  # reorder coef matrix
  m <- reorder(m, reorder.coef)

  # create output table with significance stars etc.
  ci <- logical()
  for (i in 1:length(models)) {
    if (length(models[[i]]@se) == 0 && length(models[[i]]@ci.up) > 0) {
      ci[i] <- TRUE
    } else {
      ci[i] <- FALSE
    }
  }

  # write coefficient rows
  if (output.type[1] == "ascii") {
    neginfstring <- "-Inf"
    posinfstring <- "Inf"
    se.prefix <- " ("
    se.suffix <- ")"
    star.prefix <- " "
    star.suffix <- ""
    dcolumn <- TRUE
    bold.prefix <- ""
    bold.suffix <- ""
    semicolon <- "; "
  } else if (output.type[1] == "latex") {
    neginfstring <- "\\multicolumn{1}{c}{$-\\infty$}"
    posinfstring <- "\\multicolumn{1}{c}{$\\infty$}"
    se.prefix <- " \\; ("
    se.suffix <- ")"
    star.prefix <- "^{"
    star.suffix <- "}"
    bold.prefix <- "\\mathbf{"
    bold.suffix <- "}"
    semicolon <- ";\\ "
  } else if (output.type[1] == "html") {
    neginfstring <- "-Inf"
    posinfstring <- "Inf"
    se.prefix <- " ("
    se.suffix <- ")"
    star.prefix <- paste0("<sup", dots$css.sup, ">")
    star.suffix <- "</sup>"
    dcolumn <- TRUE
    bold.prefix <- "<b>"
    bold.suffix <- "</b>"
    semicolon <- "; "
  }
  if (single.row == TRUE) {
    output.matrix <- matrix(ncol = (ncol(m) / 3) + 1, nrow = nrow(m))

    # row labels
    for (i in 1:length(rownames(m))) {
      if (output.type[1] == "latex") {
        output.matrix[i, 1] <- names2latex(rownames(m)[i])
      } else {
        output.matrix[i, 1] <- rownames(m)[i]
      }
    }

    # replace R syntax
    for (i in 1:nrow(m)) {
      if (grepl("I\\(", rownames(m)[i]) == TRUE) {
        output.matrix[i, 1] <- gsub("(.*)(?:I\\()(.+)(?:\\))(.*)", "\\1\\2\\3", output.matrix[i, 1])
      }
    }

    # coefficients and standard errors
    for (i in 1:nrow(m)) { # go through rows
      j <- 1 # column in the original, merged coef table
      k <- 2 # second column of output.matrix, i.e., coefficients
      while (j <= ncol(m)) {
        if (is.na(m[i, j])) {
          output.matrix[i, k] <- ""
        } else if (m[i, j] == -Inf) {
          output.matrix[i, k] <- neginfstring
        } else if (m[i, j] == Inf) {
          output.matrix[i, k] <- posinfstring
        } else {
          se.prefix.current <- se.prefix
          se.suffix.current <- se.suffix
          if (ci[k - 1] == TRUE) { # in case of CIs, replace parentheses by square brackets
            se.prefix.current <- gsub("\\(", "[", se.prefix.current)
            se.suffix.current <- gsub("\\)", "]", se.suffix.current)
          }
          if (is.na(m[i, j + 1])) {
            se.prefix.current <- ""
            se.suffix.current <- ""
          }
          if (ci[k - 1] == FALSE) {
            std <- paste0(se.prefix.current,
                          coeftostring(m[i, j + 1], leading.zero, digits = digits),
                          se.suffix.current)
          } else {
            std <- paste0(se.prefix.current,
                          coeftostring(m[i, j + 1], leading.zero, digits = digits),
                          semicolon,
                          coeftostring(m[i, j + 2], leading.zero, digits = digits),
                          se.suffix.current)
          }

          if (ci[k - 1] == FALSE) {
            p <- get_stars(pval = m[i, j + 2],
                           stars = stars,
                           star.symbol = star.symbol,
                           symbol = symbol,
                           star.prefix = star.prefix,
                           star.suffix = star.suffix,
                           ci = ci
            )$coefficients
          } else { # significance from confidence interval
            if (is.numeric(ci.test) && !is.na(ci.test) && bold == 0 &&
                (m[i, j + 1] > ci.test || m[i, j + 2] < ci.test)) {
              p <- paste0(star.prefix, star.symbol, star.suffix)
            } else {
              p <- ""
            }
          }

          if (isTRUE(dcolumn)) {
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
          entry <- paste0(dollar,
                          bold.pref,
                          coeftostring(m[i, j], leading.zero, digits = digits),
                          bold.suff,
                          std,
                          p,
                          dollar)
          output.matrix[i, k] <- entry
        }
        k <- k + 1
        j <- j + 3
      }
    }
  } else {
    output.matrix <- matrix(ncol = (ncol(m) / 3) + 1, nrow = 2 * nrow(m))

    # row labels
    for (i in 1:length(rownames(m))) {
      if (output.type[1] == "latex"){
        output.matrix[(i * 2) - 1, 1] <- names2latex(rownames(m)[i])
      } else {
        output.matrix[(i * 2) - 1, 1] <- rownames(m)[i]
      }
      output.matrix[(i * 2), 1] <- ""
    }

    for (i in 1:nrow(m)) {
      if (grepl("I\\(", rownames(m)[i]) == TRUE) {
        output.matrix[(i * 2) - 1, 1] <- gsub("(.*)(?:I\\()(.+)(?:\\))(.*)", "\\1\\2\\3", output.matrix[(i * 2) - 1, 1])
      }
    }

    # coefficients and standard errors
    for (i in 1:nrow(m)) {  # i = row
      j <- 1  # j = column within model (from 1 to 3)
      k <- 2  # k = column in output matrix (= model number + 1)
      while (j <= ncol(m)) {
        if (is.na(m[i, j]) || is.nan(m[i, j])) {
          output.matrix[(i * 2) - 1, k] <- "" # upper coefficient row
          output.matrix[(i * 2), k] <- "" # lower se row
        } else if (m[i, j] == -Inf) {
          output.matrix[(i * 2) - 1, k] <- neginfstring # upper row
          output.matrix[(i * 2), k] <- "" # lower se row
        } else if (m[i, j] == Inf) {
          output.matrix[(i * 2) - 1, k] <- posinfstring # upper row
          output.matrix[(i * 2), k] <- "" # lower se row
        } else {
          se.prefix.current <- "("
          se.suffix.current <- ")"
          if (ci[k - 1] == TRUE) { # in case of CIs, replace parentheses by square brackets
            se.prefix.current <- "["
            se.suffix.current <- "]"
          }
          if (is.na(m[i, j + 1])) {
            se.prefix.current <- ""
            se.suffix.current <- ""
          }
          if (ci[k - 1] == FALSE) {
            p <- get_stars(pval = m[i, j + 2],
                           stars = stars,
                           star.symbol = star.symbol,
                           symbol = symbol,
                           star.prefix = star.prefix,
                           star.suffix = star.suffix
            )$coefficients
          } else { # significance from confidence interval
            if (is.numeric(ci.test) && !is.na(ci.test) && bold == 0 &&
                (m[i, j + 1] > ci.test || m[i, j + 2] < ci.test)) {
              p <- paste0(star.prefix, star.symbol, star.suffix)
            } else {
              p <- ""
            }
          }

          if (isTRUE(dcolumn)) {
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
          output.matrix[(i * 2) - 1, k] <- paste0(dollar,
                                                  bold.pref,
                                                  coeftostring(m[i, j], leading.zero, digits = digits),
                                                  bold.suff,
                                                  p,
                                                  dollar)
          if (ci[k - 1] == FALSE) {
            output.matrix[(i * 2), k] <- paste0(dollar,
                                                se.prefix.current,
                                                coeftostring(m[i, j + 1], leading.zero, digits = digits),
                                                se.suffix.current,
                                                dollar)
          } else {
            output.matrix[(i * 2), k] <- paste0(dollar,
                                                se.prefix.current,
                                                coeftostring(m[i, j + 1], leading.zero, digits = digits),
                                                semicolon,
                                                coeftostring(m[i, j + 2], leading.zero, digits = digits),
                                                se.suffix.current,
                                                dollar)
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

  # add groups to the output matrix using 'groups' argument
  if (!is.null(groups)) {
    indentation <- "    "
    prefix <- ""
    suffix <- ""
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
      if (groups[[i]][length(groups[[i]])] - length(groups[[i]]) + 1 != groups[[i]][1]) {
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
    if (groups[[length(groups)]][length(groups[[length(groups)]])] > nrow(output.matrix)) {
      stop("'groups' argument contains indices outside the table dimensions.")
    }
    for (i in length(names(groups)):1) {
      if (output.type[1] == "latex") {
        label <- paste0(prefix, names2latex(names(groups)[i]), suffix)
      } else {
        label <- paste0(prefix, names(groups)[i], suffix)
      }
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

  # save coefficient names for matrix attributes later
  coef.names <- output.matrix[output.matrix[, 1] != "", 1]

  # combine the coefficient and gof matrices vertically
  output.matrix <- rbind(output.matrix, gof.matrix)

  # reformat output matrix columns by aligning by decimal point and adding spaces
  for (j in 2:ncol(output.matrix)) {
    x <- output.matrix[, j]
    # create first.dot: max length before first dot; and create paren.length: max
    # length of parentheses, including opening and closing parentheses
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
      # fill with spaces at the beginning
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
      spaces <- paste(rep(" ", difference), collapse = "")
      x[i] <- paste(spaces, x[i], sep = "")

      # adjust indentation for SEs
      if (single.row == TRUE) {
        paren <- attributes(parentheses)$match.length[i]
        if (paren < 0) {
          paren <- 0
        }
        difference <- paren.length - paren + 1  # +1 because strsplit takes 1 away
        spaces <- paste(rep(" ", difference), collapse = "")
        components <- strsplit(x[i], " \\(")[[1]]
        if (length(components) == 2) {
          x[i] <- paste(components[1], spaces, "(", components[2], sep = "")
        }
      }
    }

    # make all CIs have equal length
    ci.lower.length <- 0
    ci.upper.length <- 0
    for (i in 1:length(x)) {
      if (grepl("\\[.+\\]", x[i])) {
        regex <- ".*\\[(.+?);[\\\"]? (.+?)\\].*"
        first <- sub(regex, "\\1", x[i])
        first <- nchar(first)
        if (first > ci.lower.length) {
          ci.lower.length <- first
        }
        last <- sub(regex, "\\2", x[i])
        last <- nchar(last)
        if (last > ci.upper.length) {
          ci.upper.length <- last
        }
      }
    }
    for (i in 1:length(x)) {
      if (grepl("\\[.+\\]", x[i])) {
        regex <- "(.*?)\\[(.+?);[\\\"]? (.+?)\\](.*?)$"

        whitespace1 <- sub(regex, "\\1", x[i])
        whitespace1 <- sub("\\s+$", "", whitespace1)
        if (nchar(whitespace1) > 0) {
          whitespace1 <- paste0(whitespace1, " ")
        }
        whitespace2 <- sub(regex, "\\4", x[i])

        first <- sub(regex, "\\2", x[i])
        difference <- ci.lower.length - nchar(first)
        zeros <- paste(rep(" ", difference), collapse = "")
        first <- paste0(zeros, first)

        last <- sub(regex, "\\3", x[i])
        difference <- ci.upper.length - nchar(last)
        zeros <- paste(rep(" ", difference), collapse = "")
        last <- paste0(zeros, last)

        x[i] <- paste0(whitespace1, "[", first, "; ", last, "]", whitespace2)
      }
    }

    # fill with spaces at the end to make them all equally long
    max.x <- max(nchar(x))
    for (i in 1:length(x)) {
      difference <- max.x - nchar(x[i])
      spaces <- paste(rep(" ", difference), collapse = "")
      x[i] <- paste0(x[i], spaces)
    }

    output.matrix[, j] <- x
  }

  # otherwise we get duplicate model names in latex and html
  if (output.type[1] == "ascii") {
    output.matrix <- rbind(c("", mod.names), output.matrix)
  }

  # add custom columns to output.matrix using the 'custom.columns' argument
  if (!is.null(custom.columns)) { # start by checking validity of arguments
    numcoef <- nrow(m)
    if (!class(custom.columns) == "list") {
      if (length(custom.columns) != numcoef) {
        stop(paste("Custom column does not match table dimensions.", numcoef, "elements expected."))
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
      stop(paste("Length of 'custom.col.pos' does not match length of 'custom.columns'."))
    }
    if (any(custom.col.pos > ncol(output.matrix) + 1)) {
      stop(paste("The table has only", ncol(output.matrix), "columns. The",
                 "'custom.col.pos' argument does not match these dimensions."))
    }
    if (0 %in% custom.col.pos) {
      stop("0 is not a valid argument in 'custom.col.pos'. The column indices start with 1.")
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
        custom.indices <- c(custom.indices,
                            rep(TRUE, length(which(custom.col.pos == i))),
                            FALSE)
      } else {
        custom.indices <- c(custom.indices, FALSE)
      }
    }
    for (i in 1:length(custom.col.pos)) {
      if ((ncol(output.matrix) + 1) == custom.col.pos[i]) {
        custom.indices <- c(custom.indices, TRUE)
      }
    }

    # combine output matrix with custom columns
    offset <- 1
    output.count <- 0
    custom.count <- 0
    temp <- matrix(character(), nrow = nrow(output.matrix), ncol = 0)
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
            newcol[j + offset, 1] <- as.character(custom.columns[[custom.count]][j])
          } else {
            newcol[(2 * j) - (1 - offset), 1] <- as.character(custom.columns[[custom.count]][j])
          }
        }
        temp <- cbind(temp, newcol)
      }
    }
    output.matrix <- temp
  }

  # trim leading and trailing white space in output matrix
  if (isTRUE(trim)) {
    output.matrix <- trimws(output.matrix)
  }

  # attributes required for printing functions
  if (isTRUE(include.attributes)) {
    attr(output.matrix, "ci") <- ci
    attr(output.matrix, "ci.test") <- ci.test
    attr(output.matrix, "gof.names") <- gof.matrix[, 1]
    attr(output.matrix, "coef.names") <- coef.names
    attr(output.matrix, "mod.names") <- mod.names
  }

  return(output.matrix)
}

# print method for texreg table strings
print.texregTable <- function(x, ...) {
  cat(x, ...)
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
                      custom.gof.rows = NULL,
                      custom.note = NULL,
                      digits = 2,
                      leading.zero = TRUE,
                      star.symbol = "*",
                      symbol = ".",
                      override.coef = 0,
                      override.se = 0,
                      override.pvalues = 0,
                      override.ci.low = 0,
                      override.ci.up = 0,
                      omit.coef = NULL,
                      reorder.coef = NULL,
                      reorder.gof = NULL,
                      ci.force = FALSE,
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
                             custom.gof.rows = custom.gof.rows,
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
                             trim = FALSE,
                             ...)

  gof.names <- attr(output.matrix, "gof.names")
  coef.names <- attr(output.matrix, "coef.names")
  mod.names <- attr(output.matrix, "mod.names")
  ci <- attr(output.matrix, "ci")
  ci.test <- attr(output.matrix, "ci.test")

  # add spaces
  for (j in 1:ncol(output.matrix)) {
    nc <- nchar(output.matrix[, j])
    width <- max(nc)
    for (i in 1:nrow(output.matrix)) {
      spaces <- paste(rep(" ", width - nc[i]), collapse = "")
      output.matrix[i, j] <- paste0(output.matrix[i, j], spaces)
    }
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
                     output = "ascii")$note

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
                   custom.gof.rows = NULL,
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
                   ci.force = FALSE,
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

  # check dcolumn vs. bold
  if (isTRUE(dcolumn) && bold > 0) {
    dcolumn <- FALSE
    msg <- paste("The dcolumn package and the 'bold' argument cannot be used at",
                 "the same time. Switching off 'dcolumn'.")
    if (length(stars) > 1 || stars == TRUE || stars != 0) {
      warning(paste(msg, "You should also consider setting stars = 0."))
    } else {
      warning(msg)
    }
  }

  # check longtable vs. sideways
  if (isTRUE(longtable) && isTRUE(sideways)) {
    sideways <- FALSE
    msg <- paste("The longtable package and sideways environment cannot be",
                 "used at the same time. You may want to use the pdflscape package.",
                 "Switching off 'sideways'.")
    warning(msg)
  }

  # check longtable vs. float.pos
  if (isTRUE(longtable) && !(float.pos %in% c("", "l", "c", "r"))) {
    float.pos <- ""
    msg <- paste("When the longtable environment is used, the 'float.pos'",
                 "argument can only take one of the 'l', 'c', 'r', or ''",
                 "(empty) values. Setting float.pos = ''.")
    warning(msg)
  }

  # check longtable vs. scalebox
  if (isTRUE(longtable) && !is.null(scalebox)) {
    scalebox <- NULL
    warning(paste("'longtable' and 'scalebox' are not compatible. Setting scalebox = NULL."))
  }

  # matrixreg produces the output matrix
  output.matrix <- matrixreg(l,
                             single.row = single.row,
                             stars = stars,
                             custom.model.names = custom.model.names,
                             custom.coef.names = custom.coef.names,
                             custom.coef.map = custom.coef.map,
                             custom.gof.names = custom.gof.names,
                             custom.gof.rows = custom.gof.rows,
                             digits = digits,
                             leading.zero = leading.zero,
                             star.symbol = "*",
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
                             output.type = "latex",
                             include.attributes = TRUE,
                             trim = TRUE,
                             ...)

  gof.names <- attr(output.matrix, "gof.names")
  coef.names <- attr(output.matrix, "coef.names")
  mod.names <- attr(output.matrix, "mod.names")
  ci <- attr(output.matrix, "ci")
  ci.test <- attr(output.matrix, "ci.test")

  lab.length <- max(nchar(c(coef.names, gof.names))) # what is the optimal width of the labels?

  # determine column types (coef or custom) and model names in the presence of custom columns
  coltypes <- customcolumnnames(mod.names, custom.columns, custom.col.pos, types = TRUE)
  mod.names <- customcolumnnames(mod.names, custom.columns, custom.col.pos, types = FALSE)

  # define columns of the table (define now, add later)
  coldef <- ""
  if (isTRUE(no.margin)) {
    margin.arg <- "@{}"
  } else {
    margin.arg <- ""
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
        dl <- compute.width(output.matrix[, i], left = TRUE,
                            single.row = single.row, bracket = separator)
        dr <- compute.width(output.matrix[, i], left = FALSE,
                            single.row = single.row, bracket = separator)
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
                     star.symbol = "*",
                     symbol = symbol,
                     ci = ci,
                     ci.test = ci.test,
                     output = "latex")$note

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

# Microsoft Word function
wordreg <- function(l,
                    file = NULL,
                    single.row = FALSE,
                    stars = c(0.001, 0.01, 0.05),
                    custom.model.names = NULL,
                    custom.coef.names = NULL,
                    custom.coef.map = NULL,
                    custom.gof.names = NULL,
                    custom.gof.rows = NULL,
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
                    ci.force = FALSE,
                    ci.force.level = 0.95,
                    ci.test = 0,
                    groups = NULL,
                    custom.columns = NULL,
                    custom.col.pos = NULL,
                    ...) {

  if (!"rmarkdown" %in% row.names(installed.packages())) {
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
                   custom.gof.rows = custom.gof.rows,
                   digits = digits,
                   leading.zero = leading.zero,
                   star.symbol = "*", # produces error if fed star.symbol itself
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
                   output.type = "ascii",
                   include.attributes = FALSE,
                   trim = FALSE,
                   ...
  )
  wd <- getwd()
  f = tempfile(fileext = ".Rmd")
  cat(file = f, "```{r, echo = FALSE}
                    knitr::kable(mat)
                    ```", append = TRUE)
  rmarkdown::render(f, output_file = paste0(wd, "/", file))
}


# Internal helpers -------------------------------------------------------------

#' Internal functions for the \pkg{texreg} package
#'
#' @name internal
NULL

# display version number and date when the package is loaded
.onAttach <- function(libname, pkgname) {
  desc  <- packageDescription(pkgname, libname)
  packageStartupMessage(
    "Version:  ", desc$Version, "\n",
    "Date:     ", desc$Date, "\n",
    "Author:   ", "Philip Leifeld (University of Glasgow)", "\n\n",
    "Please cite the JSS article in your publications ",
    '-- see citation("texreg").'
  )
}

# extract coefficients using the broom package
broom_coefficients <- function(x) {
  out <- broom::tidy(x)
  out <- out[, c("term", "estimate", "std.error", "p.value")]
  return(out)
}

# extract gof using the broom package
broom_gof <- function(x) {
  # extract
  out <- broom::glance(x)[1, ]
  gof.decimal <- sapply(out, function(k) class(k)[1]) # type inference
  gof.decimal <- ifelse(gof.decimal %in% c("integer", "logical"), FALSE, TRUE)
  out <- data.frame("gof.names" = colnames(out),
                    "gof" = as.numeric(out),
                    "gof.decimal" = gof.decimal,
                    stringsAsFactors = FALSE)
  # rename
  gof_dict <- c(
    "adj.r.squared" = "Adj.\ R$^2$",
    "deviance" = "Deviance",
    "df" = "DF",
    "df.residual" = "DF Resid.",
    "finTol" = "Tolerance",
    "isConv" = "Convergence",
    "logLik" = "Log Likelihood",
    "null.deviance" = "Deviance (Null)",
    "p.value" = "P Value",
    "r.squared" = "R$^2$",
    "sigma" = "Sigma",
    "statistic" = "Statistic"
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
customcolumnnames <- function(modelnames, custom.columns, custom.col.pos, types = FALSE) {

  # adjust arguments
  modelnames <- c("", modelnames)
  if (is.null(custom.columns)) {
    if (isFALSE(types)) {
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
      if (i == 1 && isTRUE(types)) {
        value <- "coefnames"
      } else {
        value <- "coef"
      }
      custom.types <- c(custom.types,
                        rep("customcol", length(which(custom.col.pos == i))),
                        value)
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

  if (isTRUE(types)) {
    return(custom.types)
  } else {
    return(temp)
  }
}

get.data <- function(l, ...) {
  # if a single model is handed over, put model inside a list
  if (!"list" %in% class(l)[1]) {
    l <- list(l)
  }

  # create list of texreg objects
  models <- NULL
  for (i in 1:length(l)) {
    model <- extract(l[[i]], ...)
    if (class(model) == "list") { # must be a nested list of models (e.g., systemfit)
      models <- append(models, model)
    } else { # normal case; one model
      models <- append(models, list(model))
    }
  }
  return(models)
}

# create the star note (legend) printed at the bottom of tables and the stars
# printed next to standard errors
get_stars <- function(pval = NULL, # test statistics;
                      # leave NULL if you only want the legend
                      stars = c(0.01, 0.05, 0.1),  # numeric vector of cut-offs
                      star.symbol = "*",  # character to repeat for first 3
                      # levels of significance
                      symbol = ".",  # character for 4th level of significance
                      star.prefix = "",
                      star.suffix = "",
                      ci = FALSE,
                      ci.test = NULL,
                      css.sup = NULL,
                      output = "ascii") {

  # sanity checks and prep
  if (!output %in% c("ascii", "latex", "html")) {
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
  if (is.null(css.sup) & (output == "html")) {
    stop("To write a star note in html, you must supply 'css.sup'.")
  }
  if (length(stars) == 0) {
    stars <- NULL
  }
  p_note_flag <- any(!ci) # at least one model doesn't print confidence interval
  ci_note_flag <- any(ci) # at least one model prints a confidence interval

  symbols <- c(paste(rep(star.symbol, 3), collapse = ""),
               paste(rep(star.symbol, 2), collapse = ""),
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
      ci_note <- ""
    }
    if (output == "ascii") {
      ci_symbol <- "*"
    } else if (output == "latex") {
      ci_symbol <- "$^*$"
    } else if (output == "html") {
      ci_symbol <- paste0("<sup", css.sup, ">*</sup>")
    }
  } else {  # ci not calculated for any model -> empty ci note
    ci_note <- ""
  }

  # combine p and ci notes
  if (p_note_flag && ci_note_flag && (p_note != "") && (ci_note != "")) {
    s_note <- paste0(p_note, " (or ", ci_note, ").")
  } else if (p_note_flag && (p_note != "")) {
    s_note <- p_note
  } else if (ci_note_flag && (ci_note != "")) {
    s_note <- paste0(ci_symbol, " ", ci_note, ".")
  } else {
    s_note <- ""
  }

  # stars for individual coefficients
  if (is.null(pval) | is.null(stars)) {
    p = ""
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

  out <- list("note" = s_note, "coefficients" = p)
  return(out)
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

# replace coefs, SEs, p-values, and/or CIs by custom values if provided
override <- function(models,
                     override.coef,
                     override.se,
                     override.pvalues,
                     override.ci.low,
                     override.ci.up) {
  # check validity of override arguments for p-values and SEs
  if (class(override.se) == "list" || length(override.se) > 1 || override.se[1] != 0) {
    if (length(override.pvalues) == 1 && class(override.pvalues) != "list" && override.pvalues[1] == 0) {
      warning("Standard errors were provided using 'override.se', but p-values were not replaced!")
    }
  }
  if (class(override.pvalues) == "list" || length(override.pvalues) > 1 || override.pvalues[1] != 0) {
    if (length(override.se) == 1 && class(override.se) != "list" && override.se[1] == 0) {
      warning("p-values were provided using 'override.pvalues', but standard errors were not replaced!")
    }
  }

  # replace coefs, SEs, p-values, and/or CIs by custom values if provided
  for (i in 1:length(models)) {
    # override coefficients
    if (class(override.coef) != "list" && length(override.coef) == 1 && override.coef == 0) {
      cf <- models[[i]]@coef
    } else if (class(override.coef) == "numeric" &&
               length(models) == 1 &&
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
      warning(paste0("Number of coefficients provided does not match number of ",
                     "terms in model ", i, ". Using default values."))
      cf <- models[[i]]@coef
    } else if (class(override.coef[[i]]) != "numeric") {
      warning("Coefficients provided for model", i, "are not numeric. Using default values.")
      cf <- models[[i]]@coef
    } else {
      cf <- override.coef[[i]]
    }
    models[[i]]@coef <- cf

    # override standard errors
    if (class(override.se) != "list" && length(override.se) == 1 && override.se == 0) {
      se <- models[[i]]@se
    } else if (class(override.se) == "numeric" &&
               length(models) == 1 &&
               length(override.se) == length(models[[i]]@se)) {
      se <- override.se
    } else if (class(override.se) != "list") {
      warning("SEs must be provided as a list. Using default SEs.")
      se <- models[[i]]@se
    } else if (length(override.se) != length(models)) {
      warning("Number of SEs provided does not match number of models. Using default SEs.")
      se <- models[[i]]@se
    } else if (length(models[[i]]@se) != length(override.se[[i]])) {
      warning(paste0("Number of SEs provided does not match number of ",
                     "coefficients in model ", i, ". Using default SEs."))
      se <- models[[i]]@se
    } else if (class(override.se[[i]]) != "numeric") {
      warning(paste("SEs provided for model", i, "are not numeric. Using default SEs."))
      se <- models[[i]]@se
    } else {
      se <- override.se[[i]]
    }
    models[[i]]@se <- se

    # override p-values
    if (class(override.pvalues) != "list" && length(override.pvalues) == 1 && override.pvalues == 0) {
      pval <- models[[i]]@pvalues
    } else if (class(override.pvalues) == "numeric" &&
               length(models) == 1 &&
               length(override.pvalues) == length(models[[i]]@pvalues)) {
      pval <- override.pvalues
    } else if (class(override.pvalues) != "list") {
      warning("p-values must be provided as a list. Using default p-values.")
      pval <- models[[i]]@pvalues
    } else if (length(override.pvalues) != length(models)) {
      warning("Number of p-values provided does not match number of models. Using default p-values.")
      pval <- models[[i]]@pvalues
    } else if (length(models[[i]]@se) != length(override.pvalues[[i]])) {
      # previous line: comparison with SE because p-values can be empty
      warning(paste0("Number of p-values provided does not match number of ",
                     "coefficients in model ", i, ". Using default p-values."))
      pval <- models[[i]]@pvalues
    } else if (class(override.pvalues[[i]]) != "numeric") {
      warning(paste("p-values provided for model", i, "are not numeric. Using default p-values."))
      pval <- models[[i]]@pvalues
    } else {
      pval <- override.pvalues[[i]]
    }
    models[[i]]@pvalues <- pval

    # override lower bound of confidence intervals
    if (is.null(override.ci.low)) {
      # do nothing
    } else if (class(override.ci.low) != "list" &&
               length(override.ci.low) == 1 &&
               override.ci.low == 0) {
      ci.low <- models[[i]]@ci.low
    } else if (class(override.ci.low) == "numeric" &&
               length(models) == 1 &&
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
      warning("Lower CIs provided for model", i, "are not numeric. Using default lower CIs.")
      ci.low <- models[[i]]@ci.low
    } else {
      ci.low <- override.ci.low[[i]]
    }
    models[[i]]@ci.low <- ci.low

    # upper bound of confidence intervals
    if (is.null(override.ci.up)) {
      # do nothing
    } else if (class(override.ci.up) != "list" &&
               length(override.ci.up) == 1 &&
               override.ci.up == 0) {
      ci.up <- models[[i]]@ci.up
    } else if (class(override.ci.up) == "numeric" &&
               length(models) == 1 &&
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


# Class definition -------------------------------------------------------------

# constructor for texreg objects
createTexreg <- function(coef.names, coef, se = numeric(0),
                         pvalues = numeric(0), ci.low = numeric(0), ci.up = numeric(0),
                         gof.names = character(0), gof = numeric(0), gof.decimal = logical(0),
                         model.name = character(0)) {
  new("texreg", coef.names = coef.names, coef = coef, se = se,
      pvalues = pvalues, ci.low = ci.low, ci.up = ci.up, gof.names = gof.names,
      gof = gof, gof.decimal = gof.decimal, model.name = model.name)
}

# a class for texreg objects
setClass(Class = "texreg",
         representation = representation(
           coef.names = "character",  # row names of the coefficient block
           coef = "numeric",          # the coefficients
           se = "numeric",            # standard errors
           pvalues = "numeric",       # p-values
           ci.low = "numeric",        # lower confidence interval
           ci.up = "numeric",         # upper confidence interval
           gof.names = "character",   # row names of the goodness-of-fit block
           gof = "numeric",           # goodness-of-fit statistics
           gof.decimal = "logical",   # number of decimal places for each GOF value
           model.name = "character"   # name of the model
         ),
         validity = function(object) {
           if (length(object@coef.names) != length(object@coef)) {
             stop("coef.names and coef must have the same length!")
           }
           if (length(object@pvalues) != 0 &&
               (length(object@coef.names) != length(object@pvalues) ||
                length(object@coef) != length(object@pvalues))
           ) {
             stop(paste("pvalues must have the same length as coef.names, coef,",
                        "and se, or it must have a length of zero."))
           }
           if ((length(object@ci.low) != length(object@ci.up)) ||
               (length(object@ci.low) != 0 && length(object@ci.low) !=
                length(object@coef))) {
             stop("CIs must have a length of zero or the same length as coef.")
           }
           if (length(object@gof.names) != length(object@gof)) {
             stop("gof.names and gof must have the same length!")
           }
           if (length(object@gof.decimal) != 0 &&
               (length(object@gof.decimal) != length(object@gof) ||
                length(object@gof.decimal) != length(object@gof.names))
           ) {
             stop(paste("gof.decimal must have the same length as gof.names and",
                        "gof, or it must have a length of zero."))
           }
           if (length(object@model.name) > 1) {
             stop("Only one model name can be provided.")
           }
           return(TRUE)
         }
)

# define show method for pretty output of texreg objects
setMethod(f = "show", signature = "texreg", definition = function(object) {
  if (length(object@model.name) == 1) {
    cat(paste("Model name:", object@model.name))
  }
  if (length(object@se) == 0 && length(object@ci.up) > 0) {
    coefBlock <- cbind(object@coef, object@ci.low, object@ci.up)
    colnames(coefBlock) <- c("coef.", "lower CI", "upper CI")
  } else if (length(object@se) == 0 && length(object@pvalues) == 0) {
    cat(paste("\nNo standard errors and p-values were defined for this",
              "texreg object.\n"))
    coefBlock <- cbind(object@coef)
    colnames(coefBlock) <- "coef."
  } else if (length(object@se) == 0) {
    cat(paste("\nNo standard errors were defined for this texreg object.\n"))
    coefBlock <- cbind(object@coef, object@pvalues)
    colnames(coefBlock) <- c("coef.", "p")
  } else if (length(object@pvalues) > 0) {
    cat("\n")
    coefBlock <- cbind(object@coef, object@se, object@pvalues)
    colnames(coefBlock) <- c("coef.", "s.e.", "p")
  } else {
    cat("\nNo p-values were defined for this texreg object.\n")
    coefBlock <- cbind(object@coef, object@se)
    colnames(coefBlock) <- c("coef.", "s.e.")
  }
  rownames(coefBlock) <- object@coef.names
  dec <- object@gof.decimal
  if (length(dec) == 0) {
    cat("No decimal places were defined for the GOF statistics.\n")
  }
  cat("\n")
  print(coefBlock)
  cat("\n")
  if (length(dec) == 0) {
    gofBlock <- matrix(object@gof, ncol = 1)
    colnames(gofBlock) <- "GOF"
  } else {
    gofBlock <- data.frame(object@gof,object@gof.decimal)
    colnames(gofBlock) <- c("GOF", "dec. places")
  }
  rownames(gofBlock) <- object@gof.names
  if (nrow(gofBlock) > 0) {
    print(gofBlock)
    cat("\n")
  } else {
    cat("No GOF block defined.\n")
  }
})