#' Conversion of \R Regression Output to LaTeX or HTML Tables
#'
#' \pkg{texreg} converts coefficients, standard errors, uncertainty measures,
#' and goodness-of-fit statistics of statistical models into LaTeX or HTML
#' tables or into nicely formatted screen output for the \R console. A list of
#' several models can be combined in a single table. The output is
#' customizable. New model types can be easily implemented. Confidence
#' intervals can be used instead of standard errors and p-values.
#'
#' @keywords package
#' @author Philip Leifeld
#' @seealso \code{\link{extract}} \code{\link{texreg}}
#'
#' @references Leifeld, Philip (2013). texreg: Conversion of Statistical Model
#'   Output in R to LaTeX and HTML Tables. Journal of Statistical Software
#'   55(8): 1-24. \doi{10.18637/jss.v055.i08}.
#'
#' @name texreg-package
NULL


# texreg functions -------------------------------------------------------------

#' Convert regression output to a HTML table
#'
#' Conversion of \R regression output to a HTML table.
#'
#' The \code{htmlreg} function creates HTML code. Tables in HTML format can
#' be saved with a ".html" extension and displayed in a web browser.
#' Alternatively, they can be saved with a ".doc" extension and opened in MS
#' Word for inclusion in office documents. \code{\link{htmlreg}} also works with
#' \pkg{knitr} and HTML or Markdown. Note that the \code{inline.css},
#' \code{doctype}, \code{html.tag}, \code{head.tag}, \code{body.tag}, and
#' \code{star.symbol} arguments must be adjusted for the different purposes (see
#' the description of the arguments).
#'
#' @param inline.css Should the CSS stylesheets be embedded directly in the code
#'   of the table (\code{inline.css = TRUE}), or should the CSS stylesheets be
#'   enclosed in the <head> tag, that is, separated from the table code
#'   (\code{inline.css = FALSE})? Having inline CSS code makes the code of the
#'   table more complex, but sometimes it may be helpful when only the table
#'   shall be printed, without the head of the HTML file (for example when the
#'   table is embedded in a \pkg{knitr} report). As a rule of thumb: use inline
#'   CSS if the table is not saved to a file.
#' @param doctype Should the first line of the HTML code contain the DOCTYPE
#'   definition? If \code{TRUE}, the HTML 4 TRANSITIONAL version is used. If
#'   \code{FALSE}, no DOCTYPE will be included. Omitting the DOCTYPE can be
#'   helpful when the \pkg{knitr} package is used to generate HTML code because
#'   \pkg{knitr} requires only the plain table, not the whole HTML document
#'   including the document type declaration. Including the DOCTYPE can be
#'   helpful when the code is saved to a file, for example as an MS Word
#'   document.
#' @param html.tag Should the table code (and possibly the <body> and <head>
#'   tags) be enclosed in an <html> tag? Suppressing this tag is recommended
#'   when \pkg{knitr} is used for dynamic HTML or Markdown report generation.
#'   Including this tag is recommended when the code is saved to a file, for
#'   example as an MS Word document.
#' @param head.tag Should the <head> tag (including CSS definitions and
#'   title/caption) be included in the HTML code? Suppressing this tag is
#'   recommended when \pkg{knitr} is used for dynamic HTML or Markdown report
#'   generation. Including this tag is recommended when the code is saved to a
#'   file, for example as an MS Word document.
#' @param body.tag Should the table code be enclosed in a <body> HTML tag?
#'   Suppressing this tag is recommended when \pkg{knitr} is used for dynamic
#'   HTML or Markdown report generation. Including this tag is recommended when
#'   the code is saved to a file, for example as an MS Word document.
#' @param indentation Characters used for indentation of the HTML code. By
#'   default, \code{indentation = ""} uses no indentation. Any number of spaces
#'   or characters can be used instead. For example, \code{indentation = " "}
#'   uses two spaces of (additional) indentation for each subelement.
#' @param margin The margin around the table in pixels. This determines how much
#'   space there is around the table. To remove all space around the table, set
#'   \code{table.margin = 0}.
#' @param padding The space on the left and right of each table cell in pixels.
#' @param color The color of the table, including text and rules or lines. This
#'   can be provided as a hex RGB value or as a color string that is valid in
#'   HTML (e.g., \code{"black"}).
#' @param outer.rules The line width at the top and bottom of the table in
#'   pixels. Can be \code{outer.rules = 0} to omit outer lines.
#' @param inner.rules The horizontal line width before and after the coefficient
#'   block of the table in pixels. Can be \code{outer.rules = 0} to omit inner
#'   lines.
#' @inheritParams texreg
#' @inheritParams matrixreg
#'
#' @author Philip Leifeld
#' @family texreg
#' @seealso \code{\link{texreg-package}} \code{\link{extract}}
#'
#' @references Leifeld, Philip (2013). texreg: Conversion of Statistical Model
#'   Output in R to LaTeX and HTML Tables. Journal of Statistical Software
#'   55(8): 1-24. \doi{10.18637/jss.v055.i08}.
#'
#' @examples
#' library("nlme")
#' model.1 <- lme(distance ~ age, data = Orthodont, random = ~ 1)
#' model.2 <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1)
#' htmlreg(list(model.1, model.2),
#'         file = "texreg.doc",
#'         inline.css = FALSE,
#'         doctype = TRUE,
#'         html.tag = TRUE,
#'         head.tag = TRUE,
#'         body.tag = TRUE)
#' unlink("texreg.doc")
#'
#' @export
htmlreg <- function(l,
                    file = NULL,
                    single.row = FALSE,
                    stars = c(0.001, 0.01, 0.05),
                    custom.header = NULL,
                    custom.model.names = NULL,
                    custom.coef.names = NULL,
                    custom.coef.map = NULL,
                    custom.gof.names = NULL,
                    custom.gof.rows = NULL,
                    custom.note = NULL,
                    digits = 2,
                    leading.zero = TRUE,
                    star.symbol = "&#42;",
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
                    doctype = FALSE,
                    html.tag = FALSE,
                    head.tag = FALSE,
                    body.tag = FALSE,
                    indentation = "",
                    margin = 10,
                    padding = 5,
                    color = "#000000",
                    outer.rules = 2,
                    inner.rules = 1,
                    ...) {

  # CSS definitions
  css_table.texreg <- paste0("margin: ", margin, "px",
                             ifelse(isTRUE(center), " auto", ""),
                             ";border-collapse: collapse;border-spacing: 0px;",
                             ifelse(isTRUE(caption.above), "", "caption-side: bottom;"),
                             "color: ", color, ";border-top: ", outer.rules, "px solid ", color, ";")
  css_table.texreg_thead_.rule <- paste0("border-bottom: ", inner.rules, "px solid ", color, ";")
  css_table.texreg_tbody_tr.rule_td <- paste0("border-top: ", inner.rules, "px solid ", color, ";")
  css_table.texreg_tbody_tr.bottomrule_td <- paste0("border-bottom: ", outer.rules, "px solid ", color, ";")
  css_table.texreg_tfoot_td <- paste0("font-size: 0.8em;")
  css_table.texreg_thtd <- paste0("padding-left: ", padding, "px;padding-right: ", padding, "px;")

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
                             trim = TRUE,
                             output.type = "html",
                             ...)

  gof.names <- attr(output.matrix, "gof.names")
  coef.names <- attr(output.matrix, "coef.names")
  mod.names <- attr(output.matrix, "mod.names")
  ci <- attr(output.matrix, "ci")
  ci.test <- attr(output.matrix, "ci.test")

  coltypes <- customcolumnnames(mod.names, custom.columns, custom.col.pos, types = TRUE)
  mod.names <- customcolumnnames(mod.names, custom.columns, custom.col.pos, types = FALSE)

  # replace empty cells by HTML space
  output.matrix <- apply(output.matrix, 1:2, function(x) {
    if (x == "") {
      return("&nbsp;")
    } else {
      return(x)
    }
  })

  # write table header
  if (single.row == TRUE) {
    numcols <- 2 * length(mod.names)
  } else {
    numcols <- length(mod.names)
  }

  if (doctype == TRUE) {
    doct <- paste0("<!DOCTYPE html>\n")
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
  tabdef <- paste0(h.ind,
                   b.ind,
                   ifelse(isTRUE(inline.css),
                          paste0("<table class=\"texreg\" style=\"", css_table.texreg, "\">\n"),
                          paste0("<table class=\"texreg\">\n")))

  # set caption
  if (is.null(caption)) {
    caption <- ""
  }
  if (length(caption) > 1) {
    caption <- caption[1]
    warning("'caption' is supposed to be a character object of length 1. Using only the first element.")
  }
  if (is.na(caption)) {
    caption <- ""
  }
  if (!is.character(caption)) {
    caption <- as.character(caption)
  }
  if (caption != "") {
    cap <- paste0(h.ind, b.ind, ind, "<caption>", caption, "</caption>\n")
  } else {
    cap <- ""
  }

  # HTML header with CSS definitions
  string <- doct
  if (html.tag == TRUE) {
    string <- paste0(string, "<html lang=\"en\">\n")
  }

  if (inline.css == TRUE) {
    css.header <- ""
  } else {
    css.header <- paste0(
      h.ind, d.ind, "<style>\n",
      h.ind, d.ind, ind, "table.texreg {", css_table.texreg, "}\n",
      h.ind, d.ind, ind, "table.texreg thead .rule {", css_table.texreg_thead_.rule, "}\n",
      h.ind, d.ind, ind, "table.texreg tbody tr.rule td {", css_table.texreg_tbody_tr.rule_td, "}\n",
      h.ind, d.ind, ind, "table.texreg tbody tr.bottomrule td {", css_table.texreg_tbody_tr.bottomrule_td, "}\n",
      h.ind, d.ind, ind, "table.texreg tfoot td {", css_table.texreg_tfoot_td, "}\n",
      h.ind, d.ind, ind, "table.texreg th, table.texreg td {", css_table.texreg_thtd, "}\n",
      h.ind, d.ind, "</style>\n"
    )
  }

  if (head.tag == TRUE) {
    string <- paste0(string,
                     h.ind, "<head>\n",
                     h.ind, d.ind, "<meta charset=\"utf-8\" />\n",
                     h.ind, d.ind, "<title>", caption, "</title>\n",
                     css.header,
                     h.ind, "</head>\n")
  }
  if (body.tag == TRUE) {
    string <- paste0(string, h.ind, "<body>\n")
  }
  string <- paste0(
    string,
    tabdef,
    cap,
    h.ind, b.ind, ind, "<thead>\n",
    h.ind, b.ind, ind, ind, "<tr>\n"
  )

  # specify multicolumn header
  if (!is.null(custom.header) && length(custom.header) > 0 && !any(is.na(custom.header))) {
    if (!"list" %in% class(custom.header) || length(custom.header) >= length(mod.names) || is.null(names(custom.header)) || !all(sapply(custom.header, is.numeric))) {
      stop("'custom.header' must be a named list of numeric vectors.")
    }
    ch <- unlist(custom.header)
    for (i in 1:length(ch)) {
      if (is.na(ch[i])) {
        stop("NA values are not permitted in 'custom.header'. Try leaving out the model indices that should not be included in the custom header.")
      }
      if (ch[i] %% 1 != 0) {
        stop("The model column indices in 'custom.header' must be provided as integer values.")
      }
      if (ch[i] < 1 || ch[i] >= length(mod.names)) {
        stop("The model column indices in 'custom.header' must be between 1 and the number of models.")
      }
      if (i > 1 && ch[i] <= ch[i - 1]) {
        stop("The model column indices in 'custom.header' must be strictly increasing.")
      }
    }
    ch <- paste0(h.ind, b.ind, ind, ind, ind, "<th>&nbsp;</th>\n") # multicolumn header labels; first column is empty
    counter <- 0 # keeps track of how many columns we have processed to the left of the current start column
    for (i in 1:length(custom.header)) {
      # check if there are gaps within multicolumn blocks and throw error
      if (length(custom.header[[i]]) != custom.header[[i]][length(custom.header[[i]])] - custom.header[[i]][1] + 1) {
        stop("Each item in 'custom.header' must have strictly consecutive column indices, without gaps.")
      }

      # find out corrected column indices (ignoring the coef label column) of the current model after taking into account custom columns
      numCoefCol <- 0
      numCustomCol <- 0
      for (j in 1:length(coltypes)) {
        if (coltypes[j] == "coef") {
          numCoefCol <- numCoefCol + 1
        } else if (coltypes[j] == "customcol") {
          numCustomCol <- numCustomCol + 1
        }
        if (numCoefCol == custom.header[[i]][1]) {
          break() # break the loop if we have reached the number of models so far
        }
      }
      startIndex <- numCoefCol + numCustomCol # corrected column index
      numCoefCol <- 0
      numCustomCol <- 0
      for (j in 1:length(coltypes)) {
        if (coltypes[j] == "coef") {
          numCoefCol <- numCoefCol + 1
        } else if (coltypes[j] == "customcol") {
          numCustomCol <- numCustomCol + 1
        }
        if (numCoefCol == custom.header[[i]][length(custom.header[[i]])]) {
          break() # break the loop if we have reached the number of models so far
        }
      }
      stopIndex <- numCoefCol + numCustomCol # corrected column index

      # check if there are gaps between multicolumn headers and fill up with empty cells
      if (i > 1) {
        emptycells <- custom.header[[i]][1] - custom.header[[i - 1]][length(custom.header[[i - 1]])] - 1
        if (emptycells > 0) {
          for (j in 1:emptycells) {
            ch <- paste0(ch, h.ind, b.ind, ind, ind, ind, "<th>&nbsp;</th>\n")
            counter <- counter + 1
          }
        }
      }

      # add empty cells also for custom text columns
      if (startIndex > counter + 1) {
        difference <- startIndex - (counter + 1)
        for (j in 1:difference) {
          ch <- paste0(ch, h.ind, b.ind, ind, ind, ind, "<th>&nbsp;</th>\n")
          counter <- counter + 1
        }
      }

      # add multicolumn cells (+ 1 is for the coefficient label column)
      th <- ifelse(isTRUE(inline.css),
                   paste0("<th style=\"", css_table.texreg_thtd, css_table.texreg_thead_.rule, "\" colspan=\"", stopIndex - startIndex + 1, "\">"),
                   paste0("<th class = \"rule\" colspan=\"", stopIndex - startIndex + 1, "\">"))
      ch <- paste0(ch, h.ind, b.ind, ind, ind, ind, th, names(custom.header)[i], "</th>\n")
      counter <- counter + stopIndex - startIndex + 1
    }
    string <- paste0(string, ch)
    string <- paste0(string, h.ind, b.ind, ind, ind, "</tr>\n")
    string <- paste0(string, h.ind, b.ind, ind, ind, "<tr>\n")
  }

  # specify model names (header row)
  th <- ifelse(isTRUE(inline.css), paste0("<th style=\"", css_table.texreg_thtd, "\">"), "<th>")
  for (i in 1:length(mod.names)) {
    string <- paste0(string, h.ind, b.ind, ind, ind, ind, th, ifelse(mod.names[i] == "", "&nbsp;", mod.names[i]), "</th>\n")
  }
  string <- paste0(string, h.ind, b.ind, ind, ind, "</tr>\n")
  string <- paste0(string, h.ind, b.ind, ind, "</thead>\n")
  string <- paste0(string, h.ind, b.ind, ind, "<tbody>\n")

  # write coefficients to string object
  td <- ifelse(isTRUE(inline.css), paste0("<td style=\"", css_table.texreg_thtd, "\">"), "<td>")
  for (i in 1:nrow(output.matrix)) {
    if (i == 1 || i == nrow(output.matrix) - length(gof.names) + 1) {
      if (isTRUE(inline.css)) {
        class_rule <- paste0(" style=\"", css_table.texreg_tbody_tr.rule_td, "\"")
      } else {
        class_rule <- " class=\"rule\""
      }
    } else if (i == nrow(output.matrix)) {
      if (isTRUE(inline.css)) {
        class_rule <- paste0(" style=\"", css_table.texreg_tbody_tr.bottomrule_td, "\"")
      } else {
        class_rule <- " class=\"bottomrule\""
      }
    } else {
      class_rule <- ""
    }
    string <- paste0(string, h.ind, b.ind, ind, ind, paste0("<tr", class_rule, ">\n"))
    for (j in 1:ncol(output.matrix)) {
      string <- paste0(string, h.ind, b.ind, ind, ind, ind, td, output.matrix[i, j], "</td>\n")
    }
    string <- paste0(string, h.ind, b.ind, ind, ind, "</tr>\n")

  }
  string <- paste0(string, h.ind, b.ind, ind, "</tbody>\n")

  # stars note
  snote <- get_stars_note(stars = stars,
                          star.symbol = star.symbol,
                          symbol = symbol,
                          ci = ci,
                          ci.test = ci.test,
                          output = "html")

  if (is.null(custom.note)) {
    note <- snote
  } else if (custom.note == "") {
    note <- ""
  } else {
    note <- custom.note
    note <- gsub("%stars", snote, note)
  }
  if (note != "") {
    string <- paste0(string, h.ind, b.ind, ind, "<tfoot>\n")
    string <- paste0(string,
                     h.ind,
                     b.ind,
                     ind,
                     ind,
                     "<tr>\n",
                     h.ind,
                     b.ind,
                     ind,
                     ind,
                     ind,
                     "<td",
                     ifelse(isTRUE(inline.css), paste0(" style=\"", css_table.texreg_tfoot_td, "\""), ""),
                     " colspan=\"",
                     length(mod.names),
                     "\">",
                     note,
                     "</td>\n",
                     h.ind,
                     b.ind,
                     ind,
                     ind,
                     "</tr>\n")
    string <- paste0(string, h.ind, b.ind, ind, "</tfoot>\n")
  }

  # write table footer
  string <- paste0(string, h.ind, b.ind, "</table>\n")
  if (body.tag == TRUE) {
    string <- paste0(string, h.ind, "</body>\n")
  }
  if (html.tag == TRUE) {
    string <- paste0(string, "</html>")
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

#' Create a huxtable object from multiple statistical models
#'
#' Create a huxtable object from multiple statistical models.
#'
#' The \code{\link{huxtablereg}} function creates a
#' \code{\link[huxtable]{huxtable}} object using the \pkg{huxtable} package.
#' This allows output to HTML, LaTeX, Word, Excel, Powerpoint, and RTF. The
#' object can be formatted using \pkg{huxtable} package functions. See also
#' \code{\link[huxtable]{huxreg}}.
#'
#' @inheritParams matrixreg
#'
#' @author David Hugh-Jones
#'
#' @family texreg
#' @seealso \code{\link{texreg-package}} \code{\link{extract}}
#'
#' @examples
#' library("nlme")
#' model.1 <- lme(distance ~ age, data = Orthodont, random = ~ 1)
#' model.2 <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1)
#' if (requireNamespace("huxtable")) {
#'   hr <- huxtablereg(list(model.1, model.2))
#'   hr <- huxtable::set_bottom_border(hr, 1, -1, 0.4)
#'   hr <- huxtable::set_bold(hr, 1:nrow(hr), 1, TRUE)
#'   hr <- huxtable::set_bold(hr, 1, -1, TRUE)
#'   hr <- huxtable::set_all_borders(hr, 4, 2, 0.4)
#'   hr <- huxtable::set_all_border_colors(hr, 4, 2, "red")
#'   hr
#'   \dontrun{
#'   huxtable::quick_pdf(hr)
#'   huxtable::quick_docx(hr)
#'   # or use in a knitr document
#'   }
#' }
#'
#' @export
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
                        star.symbol = "*",
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

  mr.call <- match.call(expand.dots = TRUE)
  mr.call[[1L]] <- quote(texreg::matrixreg)
  mr.call$include.attributes <- TRUE
  mr.call$trim <- TRUE
  mx <- eval.parent(mr.call)

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

#' Flexibly choose the right table output format for use with \pkg{knitr}
#'
#' Flexibly choose the right table output format for use with \pkg{knitr}.
#'
#' This function automatically selects the right function (\link{texreg},
#' \link{screenreg}, \link{htmlreg}, or \link{matrixreg}) with the right set of
#' arguments for use with the \pkg{knitr} package, for example in RStudio. The
#' advantage of using this function with \pkg{knitr} is that the user does not
#' need to replace the \link{texreg}, \link{htmlreg} etc. function call in the
#' document when a different output format is selected.
#'
#' \link{knitreg} works with...
#' \itemize{
#'   \item \R HTML documents (\code{.Rhtml} extension)
#'   \item \R Sweave documents (\code{.Rnw} extension) for PDF output via LaTeX,
#'     rendered using...
#'     \itemize{
#'       \item the \pkg{knitr} package
#'       \item the \pkg{Sweave} package
#'     }
#'   \item \R Markdown documents (\code{.Rmd} extension), rendered as...
#'     \itemize{
#'       \item HTML documents
#'       \item PDF documents
#'       \item Word documents
#'       \item Powerpoint presentations
#'       \item Presentations (\code{.Rpres} extension, not \code{.Rmd})
#'     }
#'   \item \R Notebooks, including preview
#' }
#'
#' If Markdown and HTML rendering are selected, \link{htmlreg} arguments
#' \code{doctype = FALSE} and \code{star.symbol = "&#42;"} are set to enable
#' compatibility with Markdown. With \R HTML documents (but not Markdown) or
#' presentations (\code{.Rpres} extension), only \code{doctype = FALSE} is set.
#'
#' For PDF/LaTeX documents, the \link{texreg} argument
#' \code{use.packages = FALSE} is set to suppress any package loading
#' instructions in the preamble. The user must load any packages manually in the
#' preamble of the document.
#'
#' The \pkg{knitr} and \pkg{rmarkdown} packages must be installed for this
#' function to work.
#'
#' @param ... Arguments to be handed over to the \link{texreg}, \link{htmlreg},
#'   \link{screenreg}, or \link{matrixreg} function. See the respective help
#'   page for details.
#' @return A table as a \code{character} string in the respective output format.
#'
#' @author Philip Leifeld, with input from David Hugh-Jones
#' @family texreg
#' @seealso \code{\link{texreg-package}} \code{\link{extract}}
#'
#' @examples
#' require("nlme")
#' model.1 <- lme(distance ~ age, data = Orthodont, random = ~ 1)
#' model.2 <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1)
#' knitreg(list(model.1, model.2), center = FALSE, caption = "", table = FALSE)
#'
#' @export
knitreg <- function(...) {
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("knitreg requires the 'knitr' package to be installed.\n",
         "To do this, enter 'install.packages(\"knitr\")'.")
  }
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("knitreg requires the 'rmarkdown' package to be installed.\n",
         "To do this, enter 'install.packages(\"rmarkdown\")'.")
  }
  of <- knitr::opts_knit$get("out.format")
  if (is.null(of)) { # R Notebook preview (rendered on the R console)
    screenreg(...)
  } else if (of == "markdown") { # R Markdown document with extension .Rmd, which can be rendered to HTML, PDF, or Word
    output <- rmarkdown::all_output_formats(knitr::current_input())[1]
    if (is.null(output)) { # .Rpres presentation (not .Rmd)
      htmlreg(..., doctype = FALSE) # do not include document type because inline table
    } else if (output %in% c("html_document", "bookdown::html_document2")) { # .Rmd with HTML rendering via the rmarkdown package
      htmlreg(..., doctype = FALSE, star.symbol = "&#42;") # the star symbol must be escaped in Markdown
    } else if (output %in% c("pdf_document", "bookdown::pdf_document2", "bookdown::pdf_book")) { # .Rmd with PDF LaTeX rendering via the rmarkdown package
      texreg(..., use.packages = FALSE) # do not print \usepackage{dcolumn} etc.
    } else if (output %in% c("word_document", "powerpoint_presentation", "bookdown::word_document2")) { # .Rmd with Word/Powerpoint rendering through the rmarkdown package
      mr <- matrixreg(..., output.type = "ascii", include.attributes = FALSE, trim = TRUE)
      colnames(mr) <- mr[1, ] # set column names because we want 'kable' to draw a horizontal  line under the model names
      mr <- mr[-1, ] # remove the first row because we already set the model names as column names
      knitr::kable(mr) # use the kable function to render the table in the Word document
    } else { # unknown other output format through the rmarkdown package
      htmlreg(..., doctype = FALSE)
    }
  } else if (of == "html") { # R HTML document with extension .Rhtml (do not escape '*' symbol!)
    htmlreg(..., doctype = FALSE)
  } else if (of == "latex") { # knitr LaTeX .Rnw documents rendered to PDF
    texreg(..., use.packages = FALSE) # do not print \usepackage{dcolumn} etc.
  } else if (of == "sweave") { # Sweave LaTeX .Rnw documents rendered to PDF
    texreg(..., use.packages = FALSE) # do not print \usepackage{dcolumn} etc.
  } else if (of == "jekyll") { # not sure how Jekyll works, but I'll assume plain ASCII for now
    screenreg(...)
  } else { # whatever else knitr is throwing our way should prompt an ASCII table
    screenreg(...)
  }
}

#' Convert regression output to a \code{character} matrix
#'
#' Conversion of \R regression output to a \code{character} matrix.
#'
#' The \code{matrixreg} function creates a \code{character} matrix with the row
#' names for the coefficients and goodness-of-fit statistics in the first
#' column. The function is used under the hood by other functions like
#' \code{\link{screenreg}} or \code{\link{texreg}} but can also be called
#' directly.
#'
#' @param l A statistical model or a list of statistical models. Lists of
#'   models can be specified as \code{l = list(model.1, model.2, ...)}.
#'   Different object types can also be mixed.
#' @param single.row By default, a model parameter takes up two lines of the
#'   table: the standard error is listed in parentheses under the coefficient.
#'   This saves a lot of horizontal space on the page and is the default table
#'   format in most academic journals. If \code{single.row = TRUE} is activated,
#'   however, both coefficient and standard error are placed in a single table
#'   cell in the same line.
#' @param stars The significance levels to be used to draw stars. Between 0 and
#'   4 threshold values can be provided as a numeric vector. For example,
#'   \code{stars = numeric(0)} will not print any stars and will not print any
#'   note about significance levels below the table. \code{stars = 0.05} will
#'   attach one single star to all coefficients where the p value is below 0.05.
#'   \code{stars = c(0.001, 0.01, 0.05, 0.1)} will print one, two, or three
#'   stars, or a symbol as specified by the \code{symbol} argument depending on
#'   the p-values.
#' @param custom.model.names A character vector of labels for the models. By
#'   default, the models are named "Model 1", "Model 2", etc. Specifying
#'   \code{model.names = c("My name 1", "My name 2")} etc. overrides the default
#'   behavior.
#' @param custom.coef.names By default, \pkg{texreg} uses the coefficient names
#'   which are stored in the models. The \code{custom.coef.names} argument can
#'   be used to replace them by other character strings in the order of
#'   appearance. For example, if a table shows a total of three different
#'   coefficients (including the intercept), the argument
#'   \code{custom.coef.names = c("Intercept", "variable 1", "variable 2")} will
#'   replace their names in this order.
#'
#'   Sometimes it happens that the same variable has a different name in
#'   different models. In this case, the user can use this function to assign
#'   identical names. If possible, the rows will then be merged into a single
#'   row unless both rows contain values in the same column.
#'
#'   Where the argument contains an \code{NA} value, the original name of the
#'   coefficient is kept. For example, \code{custom.coef.names = c(NA, "age",
#'   NA)} will only replace the second coefficient name and leave the first and
#'   third name as they are in the original model.
#'
#'   See also \code{custom.coef.map} for an easier and more comprehensive way to
#'   rename, omit, and reorder coefficients.
#' @param custom.coef.map The \code{custom.coef.map} argument can be used to
#'   select, omit, rename, and reorder coefficients.
#'
#'   Users must supply a named list of this form: \code{list("x" = "First
#'   variable", "y" = NA, "z" = "Third variable")}. With that particular example
#'   of \code{custom.coef.map},
#'   \enumerate{
#'    \item coefficients will be presented in order: \code{"x"}, \code{"y"},
#'      \code{"z"}.
#'    \item variable \code{"x"} will appear as \code{"First variable"}, variable
#'      \code{"y"} will appear as \code{"y"}, and variable \code{"z"} will
#'      appear as "Third variable".
#'    \item all variables not named \code{"x"}, \code{"y"}, or \code{"z"} will
#'      be omitted from the table.
#'   }
#' @param custom.gof.names A character vector which is used to replace the
#'   names of the goodness-of-fit statistics at the bottom of the table. The
#'   vector must have the same length as the number of GOF statistics in the
#'   final table. The argument works like the \code{custom.coef.names} argument,
#'   but for the GOF values. \code{NA} values can be included where the original
#'   GOF name should be kept.
#' @param custom.gof.rows A named list of vectors for new lines at the
#'   beginning of the GOF block of the table. For example, \code{list("Random
#'   effects" = c("YES", "YES", "NO"), Observations = c(25, 25, 26))} would
#'   insert two new rows into the table, at the beginning of the GOF block
#'   (i.e., after the coefficients). The rows can contain integer, numeric, or
#'   \code{character} objects. Note that this argument is processed after the
#'   \code{custom.gof.names} argument (meaning \code{custom.gof.names} should
#'   not include any of the new GOF rows) and before the \code{reorder.gof}
#'   argument (meaning that the new GOF order specified there should contain
#'   values for the new custom GOF rows). Arguments for custom columns are not
#'   affected because they only insert columns into the coefficient block.
#' @param digits Set the number of decimal places for coefficients, standard
#'   errors and goodness-of-fit statistics. Do not use negative values! The
#'   argument works like the \code{digits} argument in the
#'   \code{\link[base:Round]{round}} function of the \pkg{base} package.
#' @param leading.zero Most journals require leading zeros of coefficients and
#'   standard errors (for example, \code{0.35}). This is also the default texreg
#'   behavior. Some journals, however, require omission of leading zeros (for
#'   example, \code{.35}). This can be achieved by setting \code{leading.zero =
#'   FALSE}.
#' @param star.symbol Alternative characters for the significance stars can be
#'   specified. This is useful if \pkg{knitr} and Markdown are used for HTML
#'   report generation. In Markdown, asterisks or stars are interpreted as
#'   special characters, so they have to be escaped. To make a HTML table
#'   compatible with Markdown, specify \code{star.symbol = "&#42;"}. Note that
#'   some other modifications are recommended for usage with \pkg{knitr} in
#'   combination with Markdown or HTML (see the \code{inline.css},
#'   \code{doctype}, \code{html.tag}, \code{head.tag}, and \code{body.tag}
#'   arguments in the \code{\link{htmlreg}} function).
#' @param symbol If four threshold values are handed over to the \code{stars}
#'   argument, p-values smaller than the largest threshold value but larger than
#'   the second-largest threshold value are denoted by this symbol. The default
#'   symbol is \code{"\\\\cdot"} for the LaTeX dot, \code{"&middot;"} for the
#'   HTML dot, or simply \code{"."} for the ASCII dot. If the
#'   \code{\link{texreg}} function is used, any other mathematical LaTeX symbol
#'   or plain text symbol can be used, for example \code{symbol = "\\\\circ"}
#'   for a small circle (note that backslashes must be escaped). If the
#'   \code{\link{htmlreg}} function is used, any other HTML character or symbol
#'   can be used. For the \code{screenreg} function, only plain text characters
#'   can be used.
#' @param override.coef Set custom values for the coefficients. New coefficients
#'   are provided as a list of numeric vectors. The list contains vectors of
#'   coefficients for each model. There must be as many vectors of coefficients
#'   as there are models. For example, if there are two models with three model
#'   terms each, the argument could be specified as \code{override.coef =
#'   list(c(0.1, 0.2, 0.3), c(0.05, 0.06, 0.07))}. If there is only one model,
#'   custom values can be provided as a plain vector (not embedded in a list).
#'   For example: \code{override.coef = c(0.05, 0.06, 0.07)}.
#' @param override.se Set custom values for the standard errors. New standard
#'   errors are provided as a list of numeric vectors. The list contains vectors
#'   of standard errors for each model. There must be as many vectors of
#'   standard errors as there are models. For example, if there are two models
#'   with three coefficients each, the argument could be specified as
#'   \code{override.se = list(c(0.1, 0.2, 0.3), c(0.05, 0.06, 0.07))}. If there
#'   is only one model, custom values can be provided as a plain vector (not
#'   embedded in a list).For example: \code{override.se = c(0.05, 0.06, 0.07)}.
#'   Overriding standard errors can be useful for the implementation of robust
#'   SEs, for example.
#' @param override.pvalues Set custom values for the p-values. New p-values are
#'   provided as a list of numeric vectors. The list contains vectors of
#'   p-values for each model. There must be as many vectors of p-values as there
#'   are models. For example, if there are two models with three coefficients
#'   each, the argument could be specified as \code{override.pvalues =
#'   list(c(0.1, 0.2, 0.3), c(0.05, 0.06, 0.07))}. If there is only one model,
#'   custom values can be provided as a plain vector (not embedded in a list).
#'   For example: \code{override.pvalues = c(0.05, 0.06, 0.07)}. Overriding
#'   p-values can be useful for the implementation of robust SEs and p-values,
#'   for example.
#' @param override.ci.low Set custom lower confidence interval bounds. This
#'   works like the other override arguments, with one exception: if confidence
#'   intervals are provided here and in the \code{override.ci.up} argument, the
#'   standard errors and p-values as well as the \code{ci.force} argument are
#'   ignored.
#' @param override.ci.up Set custom upper confidence interval bounds. This
#'   works like the other override arguments, with one exception: if confidence
#'   intervals are provided here and in the \code{override.ci.low} argument, the
#'   standard errors and p values as well as the \code{ci.force} argument are
#'   ignored.
#' @param omit.coef A character string which is used as a regular expression to
#'   remove coefficient rows from the table. For example, \code{omit.coef =
#'   "group"} deletes all coefficient rows from the table where the name of the
#'   coefficient contains the character sequence \code{"group"}. More complex
#'   regular expressions can be used to filter out several kinds of model terms,
#'   for example \code{omit.coef = "(thresh)|(ranef)"} to remove all model terms
#'   matching either \code{"thresh"} or \code{"ranef"}. The \code{omit.coef}
#'   argument is processed after the \code{custom.coef.names} argument, so the
#'   regular expression should refer to the custom coefficient names. To omit
#'   GOF entries instead of coefficient entries, use the custom arguments of the
#'   extract functions instead (see the help entry of the \code{\link{extract}}
#'   function.
#' @param reorder.coef Reorder the rows of the coefficient block of the
#'   resulting table in a custom way. The argument takes a vector of the same
#'   length as the number of coefficients. For example, if there are three
#'   coefficients, \code{reorder.coef = c(3, 2, 1)} will put the third
#'   coefficient in the first row and the first coefficient in the third row.
#'   Reordering can be sensible because interaction effects are often added to
#'   the end of the model output although they were specified earlier in the
#'   model formula. Note: Reordering takes place after processing custom
#'   coefficient names and after omitting coefficients, so the
#'   \code{custom.coef.names} and \code{omit.coef} arguments should follow the
#'   original order.
#' @param reorder.gof Reorder the rows of the goodness-of-fit block of the
#'   resulting table in a custom way. The argument takes a vector of the same
#'   length as the number of GOF statistics. For example, if there are three
#'   goodness-of-fit rows, \code{reorder.gof = c(3, 2, 1)} will exchange the
#'   first and the third row. Note: Reordering takes place after processing
#'   custom GOF names and after adding new custom GOF rows, so the
#'   \code{custom.gof.names} and \code{custom.gof.rows} arguments should follow
#'   the original order, and the \code{reorder.gof} argument should contain
#'   values for any rows that are added through the \code{custom.gof.rows}
#'   argument.
#' @param ci.force Should confidence intervals be used instead of the default
#'   standard errors and p-values? Most models implemented in the \pkg{texreg}
#'   package report standard errors and p-values by default while few models
#'   report confidence intervals. However, the functions in the \pkg{texreg}
#'   package can convert standard errors and into confidence intervals using
#'   z-scores if desired. To enforce confidence intervals instead of standard
#'   errors, the \code{ci.force} argument accepts either a logical value
#'   indicating whether all models or none of the models should be forced to
#'   report confidence intervals (\code{ci.force = TRUE} for all and
#'   \code{ci.force = FALSE} for none) or a vector of logical values indicating
#'   for each model separately whether the model should be forced to report
#'   confidence intervals (e.g., \code{ci.force = c(FALSE, TRUE, FALSE)}).
#'   Confidence intervals are computed using the standard normal distribution
#'   (z-values based on the \code{\link[stats:Normal]{qnorm}} function). The
#'   t-distribution is currently not supported because this would require each
#'   \code{\link{extract}} method to have an additional argument for the degrees
#'   of freedom.
#' @param ci.force.level If the \code{ci.force} argument is used to convert
#'   standard errors to confidence intervals, what confidence level should be
#'   used? By default, \code{0.95} is used (i.e., an alpha value of 0.05).
#' @param ci.test If confidence intervals are reported, the \code{ci.test}
#'   argument specifies the reference value to establish whether a
#'   coefficient/CI is significant. The default value \code{ci.test = 0}, for
#'   example, will attach a significance star to coefficients if the confidence
#'   interval does not contain \code{0}. A value of \code{ci.test = 1} could be
#'   useful if coefficients are provided on the odds-ratio scale, for example.
#'   If no star should be printed at all, \code{ci.test = NA} can be used. It is
#'   possible to provide a single value for all models or a vector with a
#'   separate value for each model. The \code{ci.test} argument works both for
#'   models with native support for confidence intervals and in cases where the
#'   \code{ci.force} argument is used.
#' @param bold The p-value threshold below which the coefficient shall be
#'   formatted in a bold font. For example, \code{bold = 0.05} will cause all
#'   coefficients that are significant at the 95\% level to be formatted in
#'   bold. Note that this is not compatible with the \code{dcolumn} or
#'   \code{siunitx} arguments in the \code{\link{texreg}} function. If both
#'   \code{bold} and \code{dcolumn} or \code{siunitx} are \code{TRUE},
#'   \code{dcolumn} and \code{siunitx} are switched off, and a warning message
#'   appears. Note also that it is advisable to use \code{stars = FALSE}
#'   together with the \code{bold} argument because having both bolded
#'   coefficients and significance stars usually does not make any sense.
#' @param groups This argument can be used to group the rows of the table into
#'   blocks. For example, there could be one block for hypotheses and another
#'   block for control variables. Each group has a heading, and the row labels
#'   within a group are indented. The partitions must be handed over as a list
#'   of named numeric vectors, where each number is a row index and each name is
#'   the heading of the group. Example: \code{groups = list("first group" = 1:4,
#'   "second group" = 7:8)}.
#' @param custom.columns An optional list of additional text columns to be
#'   inserted into the coefficient block of the table, for example coefficient
#'   types. The list should contain one or more character vectors with as many
#'   character or numeric elements as there are coefficients/model terms. If the
#'   vectors in the list are named, the names are used as labels in the table
#'   header. For example,
#'   \code{custom.columns = list(type = c("a", "b", "c"), 1:3)} will add two
#'   columns; the first one is labeled while the second one is not. Note that
#'   the numeric elements of the second column will be converted to
#'   \code{character} objects in this example. The consequence is that decimal
#'   alignment with the \pkg{dcolumn} package is switched off in these columns.
#'   Note that this argument is processed after any arguments that affect the
#'   number of rows.
#' @param custom.col.pos An optional integer vector of positions for the columns
#'   given in the \code{custom.columns} argument. For example, if there are
#'   three custom columns, \code{custom.col.pos = c(1, 3, 3)} will insert the
#'   first custom column before the first column of the original table and the
#'   remaining two custom columns after the second column of the original table.
#'   By default, all custom columns are placed after the first column, which
#'   usually contains the coefficient names.
#' @param dcolumn Use the \pkg{dcolumn} LaTeX package to get a nice alignment of
#'   the coefficients at the decimal separator (recommended for use with the
#'   \code{\link{texreg}} function). Note that only one of the three arguments
#'   \code{bold}, \code{dcolumn}, and \code{siunitx} can be used at a time as
#'   they are mutually incompatible.
#' @param siunitx Use the \pkg{siunitx} LaTeX package to get a nice alignment of
#'   the coefficients at the decimal separator (recommended for use with the
#'   \code{\link{texreg}} function). Note that only one of the three arguments
#'   \code{bold}, \code{dcolumn}, and \code{siunitx} can be used at a time as
#'   they are mutually incompatible.
#' @param output.type Which type of output should be produced? Valid values are
#'   \code{"ascii"} (for plain text tables), \code{"latex"} (for LaTeX markup)
#'   in the resulting table), and \code{"html"} (for HTML markup in the
#'   resulting table).
#' @param include.attributes Add some attributes to the return object for
#'   confidence intervals, coefficient names, GOF statistic names, and model
#'   names? These are used by \code{\link{texreg}} and other functions for table
#'   construction.
#' @param trim Trim leading and trailing white space in the table cells? If
#'   \code{FALSE}, the values in each column will be aligned at the decimal
#'   point, and spaces are used to make all cells equally long. This is useful
#'   for on-screen output.
#' @param ... Custom options to be passed on to the \code{\link{extract}}
#'   function. For example, most extract methods provide custom options for the
#'   inclusion or exclusion of specific goodness-of-fit statistics. See the help
#'   entries of \code{\link{extract}} for more information.
#' @return A \code{character} matrix with the coefficients and goodness-of-fit
#'   statistics and their column names.
#'
#' @author Philip Leifeld
#' @family texreg
#' @seealso \code{\link{texreg-package}} \code{\link{extract}}
#'   \code{\link{texreg}}
#'
#' @import stats
#' @export
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
                      siunitx = FALSE,
                      output.type = c("ascii", "latex", "html"),
                      include.attributes = FALSE,
                      trim = FALSE,
                      ...) {

  # check stars
  if (is.null(stars) || length(stars) == 0 || (length(stars) == 1 && is.na(stars))) {
    stars <- 0
  }
  if (any(is.na(stars)) || !is.numeric(stars) || length(stars) > 4) {
    stop("'stars' must be numeric and of length between 1 and 4.")
  }

  # unnamed arguments to environment
  dots <- list(...)

  # extract coefficients, SEs, GOFs etc. from the models and save in a list
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
        r <- paste0("<sup>2</sup>")
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
      z <- stats::qnorm(1 - ((1 - ci.force.level) / 2))
      upper <- models[[i]]@coef + (z * models[[i]]@se)
      lower <- models[[i]]@coef - (z * models[[i]]@se)
      models[[i]]@ci.low <- lower
      models[[i]]@ci.up <- upper
      models[[i]]@se <- numeric(0)
      models[[i]]@pvalues <- numeric(0)
    }
  }

  # check (and adjust) ci.test argument
  if (is.null(ci.test) || (length(ci.test) == 1 && is.na(ci.test))) {
    ci.test <- NA
  } else if (!is.numeric(ci.test) && !all(is.na(ci.test))) {
    stop("'ci.test' must be numeric or NA.")
  }
  if (length(ci.test) != 1 && length(ci.test) != length(models)) {
    stop("'ci.test' must be either a single value for all models or one value per model or NA.")
  }
  if (length(ci.test) == 1) {
    ci.test <- rep(ci.test, length(models))
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
  if (dcolumn == FALSE && siunitx == FALSE && output.type[1] == "latex") {
    dollar <- "$"
  } else {
    dollar <- ""
  }
  gof.matrix <- matrix(nrow = nrow(gofs), ncol = ncol(gofs) + 1)  # including labels
  if (nrow(gof.matrix) > 0) {
    for (i in 1:nrow(gofs)) {
      gof.matrix[i, 1] <- rownames(gofs)[i]
      for (j in 1:ncol(gofs)) {
        strg <- coeftostring(gofs[i, j], leading.zero, digits = decimal.matrix[i, j])
        gof.matrix[i, j + 1] <- paste0(dollar, strg, dollar)
      }
    }
  }

  # add custom GOF rows
  if (!is.null(custom.gof.rows)) {
    if (!"list" %in% class(custom.gof.rows)) {
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
        if (!is.numeric(custom.gof.rows[[i]])) { # put NA for character objects, for example fixed effects
          dec <- NA
        } else if (all(custom.gof.rows[[i]] %% 1 == 0)) { # put 0 for integers
          dec <- 0
        } else { # put the respective decimal places if numeric but not integer
          dec <- digits
        }
        newValues <- sapply(custom.gof.rows[[i]], function(x) { # format the different values of the new row
          if (is.character(x) && output.type[1] == "latex" && isTRUE(dcolumn)) {
            paste0("\\multicolumn{1}{c}{", x, "}")
          } else if (is.character(x) && output.type[1] == "latex" && isTRUE(siunitx)) {
            paste0("{", x, "}")
          } else if (is.character(x) && output.type[1] == "latex" && !isTRUE(dcolumn) && !isTRUE(siunitx)) {
            coeftostring(x, leading.zero, digits = dec) # omit dollars around value if character and no dcolumn/siunitx (would otherwise be printed in italics)
          } else {
            paste0(dollar, coeftostring(x, leading.zero, digits = dec), dollar)
          }
        })
        gof.matrix <- rbind(c(names(custom.gof.rows)[i], newValues), gof.matrix) # insert GOF name and GOF values as new row
      }
    }
  }
  colnames(gof.matrix) <- NULL
  if (output.type[1] == "latex") { # if LaTeX, replace underscores etc in GOF labels
    gof.matrix[, 1] <- names2latex(gof.matrix[, 1])
  }

  # apply custom coefficient map using 'custom.coef.map' argument
  if (!is.null(custom.coef.map)) {
    # sanity checks
    if (!"list" %in% class(custom.coef.map) || is.null(names(custom.coef.map))) {
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
    m <- out
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
    custom.coef.names[is.na(custom.coef.names)] <- rownames(m)[is.na(custom.coef.names)]
    m <- m[idx, , drop = FALSE]
    row.names(m) <- custom.coef.names
  }

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
    rows <- matrix(NA, nrow = 0, ncol = orig.width) # create matrix where re-arranged rows will be stored
    for (j in 1:orig.width) {                       # go through columns in m
      current.name <- unique.names[i]               # save row name
      nonNa <- m[rownames(m) == current.name, j]    # create a vector of values with same rowname in the col
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
    } else if (!is.null(names(l)) && !is.na(names(l)[i]) && names(l)[i] != "" && "list" %in% class(l)) {
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

  # process star symbols; create 'symbols' vector
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
    star.prefix <- paste0("<sup>")
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

          if (ci[k - 1] == FALSE) { # attach SE stars to the matrix cell
            pval <- m[i, j + 2]
            if (is.null(pval) || is.null(stars) || length(stars) == 0) {
              p = ""
            } else {
              if (is.na(pval)) {
                pval <- 1.0
              }
              st <- sort(stars)
              idx <- (pval < st)
              if (any(idx)) {  # choose lowest threshold (relies on previous sorting)
                p <- paste0(star.prefix, symbols[idx][1], star.suffix)
              } else {  # not significant
                p <- ""
              }
            }
          } else { # significance from confidence interval
            if (!is.na(ci.test[k - 1]) && bold == 0 &&
                (m[i, j + 1] > ci.test[k - 1] || m[i, j + 2] < ci.test[k - 1])) {
              p <- paste0(star.prefix, star.symbol, star.suffix)
            } else {
              p <- ""
            }
          }

          if (isTRUE(dcolumn) || isTRUE(siunitx)) {
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
          if (ci[k - 1] == FALSE) { # attach SE stars to the matrix cell
            pval <- m[i, j + 2]
            if (is.null(pval) || is.null(stars) || length(stars) == 0) {
              p = ""
            } else {
              if (is.na(pval)) {
                pval <- 1.0
              }
              st <- sort(stars)
              idx <- (pval < st)
              if (any(idx)) {  # choose lowest threshold (relies on previous sorting)
                p <- paste0(star.prefix, symbols[idx][1], star.suffix)
              } else {  # not significant
                p <- ""
              }
            }
          } else { # significance from confidence interval
            if (!is.na(ci.test[k - 1]) && bold == 0 &&
                (m[i, j + 1] > ci.test[k - 1] || m[i, j + 2] < ci.test[k - 1])) {
              p <- paste0(star.prefix, star.symbol, star.suffix)
            } else {
              p <- ""
            }
          }

          if (isTRUE(dcolumn) || isTRUE(siunitx)) {
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
      output.matrix <- output.matrix[-se.missing, , drop = FALSE]
    }
  }

  # add groups to the output matrix using 'groups' argument
  if (!is.null(groups)) {
    if (output.type[1] == "latex") {
      indentation <- "\\quad "
    } else if (output.type[1] == "html") {
      indentation <- "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"
    } else {
      indentation <- "    "
    }
    prefix <- ""
    suffix <- ""
    if (!"list" %in% class(groups)) {
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
    if (!"list" %in% class(custom.columns)) {
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
        for (j in 1:numcoef) {
          if (single.row == TRUE) {
            j_index <- j
          } else {
            j_index <- (2 * j) - 1
          }
          if (output.type[1] == "ascii") { # ASCII tables have the model names already built in, so the row counter needs to be corrected
            j_index <- j_index + 1
          }
          newcol[j_index, 1] <- as.character(custom.columns[[custom.count]][j])
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

  # replace model names for ASCII tables because they already have the model names but fail to include custom column names
  if (output.type[1] == "ascii") {
    mod.names <- customcolumnnames(mod.names, custom.columns, custom.col.pos, types = FALSE)
    output.matrix[1, ] <- mod.names
  }

  return(output.matrix)
}

#' Prints a \code{texregTable} object.
#'
#' @param x A \code{texregTable} argument, as produced by \code{\link{texreg}}
#'   and related functions.
#' @param ... Additional arguments for the \code{\link[base]{cat}} function.
#'
#' @method print texregTable
#'
#' @author Philip Leifeld
#'
#' @export
print.texregTable <- function(x, ...) {
  cat(x, ...)
}

#' Create coefficient plots from statistical model output using \pkg{ggplot2}.
#'
#' Create coefficient plots of \R regression output using \pkg{ggplot2}.
#'
#' The \code{plotreg} function produces coefficient plots (i.e., forest plots
#' applied to point estimates and confidence intervals) and works much like the
#' \code{\link{screenreg}}, \code{\link{texreg}}, \code{\link{htmlreg}},
#' \code{\link{matrixreg}} and \code{\link{wordreg}} functions. It accepts a
#' single model or multiple statistical models as input and internally extracts
#' the relevant data from the models. If confidence intervals are not defined in
#' the extract method of a statistical model (see \link{extract}), the default
#' standard errors are converted to confidence intervals. Most of the arguments
#' work like in the \code{\link{screenreg}}, \code{\link{texreg}}, and
#' \code{\link{htmlreg}} \code{\link{matrixreg}}, and \code{\link{wordreg}}
#' functions. It is possible to display the plots in two ways: using the
#' \code{type = "facet"} argument, one forest plot applied to point estimates
#' and confidence intervals will be visualized in case there is only one model.
#' If there is more than one model, each one will be plotted next to the other
#' as a separate facet; using the \code{type = "forest"} argument, coefficients
#' from one or more models will be grouped together and displayed as a single
#' forest plot.
#'
#' @param custom.title With this argument, a replacement text for the
#'   \code{ggtitle}, which provides a title above the diagram, can be provided.
#'   If an empty character object is provided (\code{custom.title = ""}), the
#'   title will be omitted completely.
#' @param override.pval Set custom values for the p-values. New p-values are
#'   provided as a list of numeric vectors. The list contains vectors of
#'   p-values for each model. There must be as many vectors of p-values as there
#'   are models. For example, if there are two models with three coefficients
#'   each, the argument could be specified as \code{override.pvalues =
#'   list(c(0.1, 0.2, 0.3), c(0.05, 0.06, 0.07))}. If there is only one model,
#'   custom values can be provided as a plain vector (not embedded in a list).
#'   For example: \code{override.pvalues = c(0.05, 0.06, 0.07)}. Overriding
#'   p-values can be useful for the implementation of robust SEs and p-values,
#'   for example.
#' @param ci.level If standard errors are converted to confidence intervals
#'   (because a model does not natively support CIs), what confidence level
#'   should be used for the outer confidence interval? By default, \code{0.95}
#'   is used (i.e., an alpha value of 0.05).
#' @param ci.test If confidence intervals are reported, the \code{ci.test}
#'   argument specifies the reference value to establish whether a
#'   coefficient/CI is significant. The default value \code{ci.test = 0}, for
#'   example, will display coefficients with a round circle and the red color
#'   if the confidence interval does not contain \code{0}. A value of
#'   \code{ci.test = 1} could be useful if coefficients are provided on the
#'   odds-ratio scale, for example. It is possible to provide a single value
#'   for all models or a vector with a separate value for each model
#'   (even if it would make the plot hard to read). The \code{ci.test} argument
#'    works both for models with native support for confidence intervals and
#'    in cases where the \code{ci.force} argument is used.
#' @param type The default option is \code{type = "facet"}. If only one model is
#'   specified, it will print one forest plot applied to point estimates and
#'   confidence intervals. If more than one model is specified, it will print
#'   as many facets as the number of models in a column of plots. Alternatively,
#'   if \code{type = "forest"} is specified, coefficients from one or more
#'   models will be grouped together and displayed as a single forest plot.
#' @param theme The \code{theme} argument can be used to customize the
#'   appearance of the plot. The default theme is \code{theme_bw}. It can be
#'   replaced by any other \pkg{ggplot2} theme. See
#'   \code{\link[ggplot2]{ggtheme}} for details.
#' @param signif.light Color of outer confidence intervals for significant model
#'   terms.
#' @param signif.medium Color of inner confidence intervals for significant
#'   model terms.
#' @param signif.dark Color of point estimates and labels for significant model
#'   terms.
#' @param insignif.light Color of outer confidence intervals for insignificant
#'   model terms.
#' @param insignif.medium Color of inner confidence intervals for insignificant
#'   model terms.
#' @param insignif.dark Color of point estimates and labels for insignificant
#'   model terms.
#'
#' @return Coefficient plot as a \pkg{ggplot2} \code{gg} object if
#'   \code{file = FALSE}. \code{NULL} otherwise.
#'
#' @inheritParams texreg
#' @inheritParams matrixreg
#'
#' @author Claudia Zucca, Philip Leifeld
#' @family texreg
#' @seealso \code{\link{texreg-package}} \code{\link{extract}}
#'   \code{\link{texreg}} \code{\link{matrixreg}}
#'
#' @examples
#' \dontrun{
#' # example from the 'lm' help file:
#' ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
#' trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
#' group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
#' weight <- c(ctl, trt)
#' lm.D9 <- lm(weight ~ group)
#' lm.D90 <- lm(weight ~ group - 1)
#' plotreg(lm.D9) # plot model output as a diagram
#'
#' # customize theme and title and save as a PDF file.
#' plotreg(lm.D9,
#'         theme = theme_dark(),
#'         ggtitle = "my title",
#'         file = "myplot.pdf")
#' unlink("myplot.pdf")
#'
#' # group coefficients from multiple models
#' plotreg(list(lm.D9, lm.D90), type = "forest")
#' }
#'
#' @import stats
#' @export
plotreg <- function(l,
                    file = NULL,
                    custom.model.names = NULL,
                    custom.title = NULL,
                    custom.coef.names = NULL,
                    custom.coef.map = NULL,
                    custom.note = NULL,
                    override.coef = 0,
                    override.se = 0,
                    override.pval = 0,
                    override.ci.low = 0,
                    override.ci.up = 0,
                    override.pvalues = 0,
                    omit.coef = NULL,
                    reorder.coef = NULL,
                    ci.level = 0.95,
                    ci.force = FALSE,
                    ci.force.level = 0.95,
                    ci.test = 0,
                    type = "facet",
                    theme = NULL,
                    signif.light = "#FBC9B9",
                    signif.medium = "#F7523A",
                    signif.dark = "#BD0017",
                    insignif.light = "#C5DBE9",
                    insignif.medium = "#5A9ECC",
                    insignif.dark = "#1C5BA6",
                    ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plotreg requires the 'ggplot2' package to be installed.\n",
         "To do this, enter 'install.packages(\"ggplot2\")'.")
  }

  if (is.null(theme)) {
    theme <- ggplot2::theme_bw()
  }

  if (!is.null(omit.coef) && !is.na(omit.coef) && !is.character(omit.coef)) {
    stop("'omit.coef' must be a character string!")
  }
  if (is.null(omit.coef)) {
    omit.coef <- NA
  }

  # Extract texreg objects and override data
  models <- get.data(l, ...)
  models <- override(models,
                     override.coef,
                     override.se,
                     override.pval,
                     override.ci.low,
                     override.ci.up)
  # Custom model names
  model.names <- character()
  if (is.null(custom.model.names)) {
    model.names <- paste("Model", 1:length(l))
  } else if (length(custom.model.names) == 1) {
    model.names <- rep(custom.model.names, length(l))
  } else if (length(custom.model.names) != length(l)) {
    stop("The 'custom.model.names' argument must have the same length as the 'l' argument.")
  } else {
    model.names <- custom.model.names
  }

  # Initialize variables to store the data for plotting
  co <- se <- co.names <- pv <- lab <- ci.low <- ci.up <- NULL

  # Custom coefficients names
  for (i in 1:length(models)) {
    if (!is.null(custom.coef.names)) {
      if ("list" %in% class(custom.coef.names)) {
        if (length(custom.coef.names[[i]]) == length(models[[i]]@coef.names)) {
          models[[i]]@coef.names <- custom.coef.names[[i]]
        } else {
          stop(paste0("Model ", i, ": wrong number of custom coefficient names."))
        }
      } else if (is.character(custom.coef.names)) {
        if (length(custom.coef.names) == length(models[[i]]@coef.names)) {
          models[[i]]@coef.names <- custom.coef.names
        } else {
          stop(paste0("Model ", i, ": wrong number of custom coefficient names."))
        }
      } else {
        stop("Custom coefficient names must be provided as a list of character vectors.")
      }
    }

    # Store data in the inizialized variables for plotting
    coef.names <- models[[i]]@coef.names
    co.names <- append(co.names, coef.names)
    coefs <- models[[i]]@coef
    co <- append(co, coefs)
    label <- rep(paste0("Model", i), length(models[[i]]@coef))
    lab <- append(lab, label)
    if (length(models[[i]]@se) > 0) {
      se <- append(se, models[[i]]@se)
    }
    if (length(models[[i]]@pvalues) > 0) {
      pv <- append(pv, models[[i]]@pvalues)
    }
    if (length(models[[i]]@ci.low) > 0) {
      ci.low <- append(ci.low, models[[i]]@ci.low)
    }
    if (length(models[[i]]@ci.up) > 0) {
      ci.up <- append(ci.up, models[[i]]@ci.up)
    }

    # Create confidence intervals using ci.force argument if it is set to TRUE
    # (code replicated from matrixreg, but modified where needed)
    if (isTRUE(ci.force[i])) {
      ci.low <- ci.up <- NULL
    }
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
      if (isTRUE(ci.force[i]) && length(models[[i]]@se) > 0) {
        z <- stats::qnorm(1 - ((1 - ci.force.level) / 2))
        upper <- models[[i]]@coef + (z * models[[i]]@se)
        lower <- models[[i]]@coef - (z * models[[i]]@se)

        # store the coefficients as a ordered vector from model 1 to N
        ci.low <- append(ci.low, lower)
        ci.up <- append(ci.up, upper)
      }
    }
  }

  # Building the main dataframe to store data for plotting
  dataframe <- data.frame(cbind(co.names, co, lab))

  if (length(se) > 0) {
    dataframe <- cbind(dataframe, se)
  }
  if (length(pv) > 0) {
    dataframe <- cbind(dataframe, pv)
  }
  if (length(ci.low) > 0) {
    dataframe <- cbind(dataframe, ci.low)
  }
  if (length(ci.up) > 0) {
    dataframe <- cbind(dataframe, ci.up)
  }

  # Impose numeric class to every variable in the dataframe
  dataframe$co <- as.numeric(as.character(dataframe$co))

  if (length(dataframe$se) > 0) {
    dataframe$se <- as.numeric(as.character(dataframe$se))
  }
  if (length(dataframe$pv) > 0) {
    dataframe$pv <- as.numeric(as.character(dataframe$pv))
  }
  if (length(dataframe$ci.low) > 0) {
    dataframe$ci.low <- as.numeric(as.character(dataframe$ci.low))
  }
  if (length(dataframe$ci.up) > 0) {
    dataframe$ci.up <- as.numeric(as.character(dataframe$ci.up))
  }

  # Remove se and pv from dataframe in case ci.up and ci.low are provided
  if (length(dataframe$ci.up > 0) && length(dataframe$ci.low > 0)) {
    dataframe <- dataframe[, !(colnames(dataframe) %in% c("se","pv"))]
  }

  # Initialize variable to store intervals for bars in plot
  # (name kept from old version where inner bars were also displayed)
  lower.outer <- upper.outer <- NULL

  # Create intervals from standard error OR assign ci.low and ci.up to the
  # variables called for plotting
  if (isFALSE(ci.force[i]) && length(dataframe$se) > 0) {
    lower.outer <- dataframe$co - 2 * dataframe$se
    upper.outer <- dataframe$co + 2 * dataframe$se
  } else if (length(dataframe$ci.low) > 0 && length(dataframe$ci.up) > 0) {
    lower.outer <- dataframe$ci.low
    upper.outer <- dataframe$ci.up
  } else if (length(dataframe$se) == 0 && length(dataframe$pvalues) > 0) {
    stop("Model has p-values but no SEs. SEs or CIs are required for plotting.")
  }

  # Bind the lower.outer and upper.outer in the dataframe
  if (length(lower.outer) > 0) {
    dataframe <- cbind(dataframe, lower.outer)
  }
  if (length(upper.outer) > 0) {
    dataframe <- cbind(dataframe, upper.outer)
  } else if (length(upper.outer) == 0 && length(lower.outer) == 0) {
    stop("Model does not have either SEs or CIs. SEs or CIs are required for plotting bars.")
  }

  # Errors for the ci.test argument (adapted from matrixreg)
  if (is.null(ci.test) || (length(ci.test) == 1 && is.na(ci.test))) {
    ci.test <- NA
  } else if (!is.numeric(ci.test) && !all(is.na(ci.test))) {
    stop("'ci.test' must be numeric or NA.")
  }
  if (length(ci.test) != 1 && length(ci.test) != length(models)) {
    stop("'ci.test' must be either a single value for all models or one value per model or NA.")
  }
  if (length(ci.test) == 1) {
    ci.test <- rep(ci.test, length(models))
  }

  # Assign significance labels to terms
  if (length(dataframe$pv) == 0 & length(dataframe$ci.low) == 0) {
    stop("Impossible to estimate significance since no CIs or pvalues are provided.")
  } else if (length(dataframe$pv) > 0) {
    # Significance from SE OR CI
    signif.outer <- dataframe$pv < (1 - ci.level)
  } else if (length(dataframe$ci.low) > 0) {
    signif.outer <- ((dataframe$ci.low > ci.test & dataframe$ci.up > ci.test) |
                       (dataframe$ci.low < ci.test & dataframe$ci.up < ci.test))
  }

  # Insert significance in dataframe
  if (is.logical(signif.outer) && length(signif.outer) == length(dataframe$co)) {
    signif <- signif.outer
    dataframe <- cbind(dataframe, signif)
  } else {
    stop("'signif.outer' does not correspond to the intervals provided.")
  }

  # Make sure lab and co.names are factors, otherwise it clashes with ggplot2
  dataframe$lab <- as.factor(dataframe$lab)
  dataframe$co.names <- as.factor(dataframe$co.names)

  # Error message for coefficient misspecification
  if (length(co) == 0) {
    stop(paste("No coefficients available. Was the 'omit.coef' argument",
               "misspecified? If coefficients were renamed using the",
               "'custom.coef.names' argument, 'omit.coef' refers to the renamed",
               "coefficients."))
  }

  # Processing omit.coef
  if (!is.na(omit.coef)) {
    dataframe <- dataframe[!grepl(omit.coef, dataframe$co.names), ]
  }

  # Processing custom coef. name
  if (length(custom.coef.names) > 0) {
    levels(dataframe$co.names) <- custom.coef.names
  }

  # Processing custom model name
  if (length(custom.model.names) > 0) {
    levels(dataframe$lab) <- custom.model.names
  }

  # Starting ggplot functions
  if (type == "facet") {

    # Processing reorder.coef
    dataframe$co.names <- ordered(dataframe$co.names, levels = (unique(as.character(dataframe$co.names))))
    if (length(reorder.coef) > 0) {
      # Reorder coef
      reorder.coef<- rev(reorder.coef)
      dataframe$co.names <- factor(dataframe$co.names, levels(dataframe$co.names)[reorder.coef])
    } else {
      # Reorder coef in original order
      startlevel<- c(length(dataframe$co.names):1)
      dataframe$co.names <- factor(dataframe$co.names, levels(dataframe$co.names)[startlevel])
    }

    # Processing custom.coef.map
    if (length(custom.coef.map) > 0) {
      # Keep only coefficients selected by the user and replace NAs if any
      idx <- is.na(custom.coef.map)
      custom.coef.map[idx] <- names(custom.coef.map)[idx]
      selectcoef <- names(custom.coef.map)
      keepcoef <- paste(selectcoef, collapse = "|")
      dataframe <- dataframe[grepl(keepcoef, dataframe$co.names), ]
      # Resize number of levels in factor to avoid error
      dataframe$co.names <- ordered(dataframe$co.names, levels = unique(as.character(dataframe$co.names)))
      # Rename selected coefficients
      renamedcoef <- paste(unlist(custom.coef.map), sep = ", ")
      levels(dataframe$co.names) <- renamedcoef
      # Invert order factors for plot
      coeforder<- c(length(dataframe$co.names):1)
      dataframe$co.names <- factor(dataframe$co.names, levels(dataframe$co.names)[coeforder])
    }

    p <- ggplot2::ggplot(dataframe, ggplot2::aes(co.names, co)) +
      ggplot2::geom_hline(yintercept = ci.test,
                          lty = 2, lwd = 1,
                          colour = "grey50") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = lower.outer, ymax = upper.outer),
                             lwd = 1,
                             colour = ifelse(sapply(dataframe$signif, isTRUE) , signif.light, insignif.light),
                             width = 0) +
      ggplot2::geom_point(size = 3,
                          pch = ifelse(sapply(dataframe$signif, isTRUE), 21, 22),
                          fill = ifelse(sapply(dataframe$signif, isTRUE), signif.dark, insignif.dark)) +
      ggplot2::coord_flip()
    p <- p + theme
    p <- p + ggplot2::xlab(" ") +
      ggplot2::facet_wrap(~ lab,
                          strip.position = "left",
                          nrow = length(dataframe),
                          scales = "free_y")

    if (length(custom.title) > 0) {
      p <- p + ggplot2::ggtitle(custom.title)
    }

  } else if (type == "forest") {
    if (length(reorder.coef) > 0) {
      # Reorder levels
      dataframe$co.names <- ordered(dataframe$co.names, levels = (unique(as.character(dataframe$co.names))))
      # Reorder coef
      dataframe$co.names <- factor(dataframe$co.names, levels(dataframe$co.names)[reorder.coef])
    } else {
      # Reorder levels
      dataframe$co.names <- ordered(dataframe$co.names, levels = (unique(as.character(dataframe$co.names))))
      # Reorder coef in original order
      startlevel<- c(1:length(dataframe$co.names))
      dataframe$co.names <- factor(dataframe$co.names, levels(dataframe$co.names)[startlevel])
    }

    if (length(custom.coef.map) > 0) {
      # Keep only coefficients selected by the user and replace NAs if any
      idx <- is.na(custom.coef.map)
      custom.coef.map[idx] <- names(custom.coef.map)[idx]
      selectcoef <- names(custom.coef.map)
      keepcoef <- paste(selectcoef, collapse = "|")
      dataframe <- dataframe[grepl(keepcoef, dataframe$co.names), ]
      # Resize number of levels in factor to avoid error
      dataframe$co.names <- ordered(dataframe$co.names, levels = unique(as.character(dataframe$co.names)))
      # Rename selected coefficients
      renamedcoef <- paste(unlist(custom.coef.map), sep = ", ")
      levels(dataframe$co.names) <- renamedcoef
      # Invert order factors for plot
      coeforder<- c(length(dataframe$lab):1)
      dataframe$lab <- factor(dataframe$lab, levels(dataframe$lab)[coeforder])
    }

    # Put models in the right order
    laborder<- c(length(dataframe$lab):1)
    dataframe$lab <- factor(dataframe$lab, levels(dataframe$lab)[laborder])

    # Plot
    p <- ggplot2::ggplot(dataframe, ggplot2::aes(lab, co)) +
      ggplot2::geom_hline(yintercept = ci.test,
                          lty = 2,
                          lwd = 1,
                          colour = "grey50") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = lower.outer, ymax = upper.outer),
                             lwd = 1,
                             colour = ifelse(sapply(dataframe$signif, isTRUE) , signif.light, insignif.light),
                             width = 0) +
      ggplot2::geom_point(size = 3,
                          pch = ifelse(sapply(dataframe$signif, isTRUE), 21, 22),
                          fill = ifelse(sapply(dataframe$signif, isTRUE), signif.dark, insignif.dark)) +
      ggplot2::coord_flip()
    p <- p + theme
    p <- p + ggplot2::xlab(" ") +
      ggplot2::facet_wrap(~ co.names,
                          strip.position = "left",
                          nrow = length(dataframe),
                          scales = "free_y")

    if (length(custom.title) > 0) {
      p <- p + ggplot2::ggtitle(custom.title)
    }
  }

  # It adds message to p as ylab; adds output meassages to paste that comes as R
  # output (not in plot)
  ci.level.perc <- ci.level*100
  if (isFALSE(ci.force[i]) && length(dataframe$se) > 0) {
    p <- p + ggplot2::ylab(paste0("Bars denote SEs (",
                                  deparse1(substitute(ci.level.perc)),
                                  "%). Circle points denote significance."))
    if (length(custom.note) > 0) {
      p <- p + ggplot2::ylab(custom.note)
    }
    if (length(models) == 1) {
      message(paste0("Model: bars denote standard errors (",
                     deparse1(substitute(ci.level.perc)), "%)."))
    } else if (length(models) > 1) {
      message(paste0("Models: bars denote standard errors (" ,
                     deparse1(substitute(ci.level.perc)), "%)."))
    }
  } else if (isTRUE(ci.force[i]) && length(dataframe$se) == 0) {
    p <- p + ggplot2::ylab(paste0("Bars denote CIs (",
                                  deparse1(substitute(ci.level.perc)),
                                  "%) computed from SEs. Circle points denote significance."))
    if (length(custom.note) > 0) {
      p <- p + ggplot2::ylab(custom.note)
    }
    if (length(models) == 1) {
      message(paste0("Model: bars denote ", ci.level,
                     " confidence intervals."))
    } else if (length(models) > 1) {
      message(paste0("Models: bars denote ", ci.level,
                     " confidence intervals."))
    }
  } else {
    p <- p + ggplot2::ylab(paste0("Bars denote CIs (",
                                  deparse1(substitute(ci.level.perc)),
                                  "%). Circle points denote significance."))
    if (length(custom.note) > 0) {
      p <- p + ggplot2::ylab(custom.note)
    }
    if (length(models) == 1) {
      message(paste0("Model: bars denote ", ci.level, " confidence intervals."))
    } else if (length(models) > 1) {
      message(paste0("Models: bars denote ", ci.level, " confidence intervals."))
    }
  }

  # Print plot in R
  if (is.null(file)) {
    return(p)
  } else {
    ggplot2::ggsave(filename = file, plot = p)
  }
}

#' Convert regression output to an ASCII table
#'
#' Conversion of \R regression output to an ASCII table for display on screen.
#'
#' The \code{\link{screenreg}} function creates text representations of tables
#' and prints them to the \R console. This is an alternative to the
#' \code{\link[base]{summary}} function and serves easy model comparison.
#' Moreover, once a table has been prepared in the \R console, it can be later
#' exported to LaTeX or HTML with little extra effort because the majority of
#' arguments of the different functions are identical.
#'
#' @param column.spacing The amount of space between any two columns of a table.
#'   By default, two spaces are used. If the tables do not fit on a single page
#'   horizontally, the value can be set to \code{1} or \code{0}.
#' @param outer.rule The character which is used to draw the outer horizontal
#'   line above and below a table. If an empty character object is provided
#'   (i.e., \code{outer.rule = ""}), there will be no outer horizontal lines.
#'   Recommended values are \code{""}, \code{"="}, \code{"-"}, \code{"_"}, or
#'   \code{"#"}.
#' @param inner.rule The character used to draw the inner horizontal line above
#'   and below a table. If an empty \code{character} object is provided (i.e.,
#'   \code{outer.rule = ""}), there will be no inner horizontal lines.
#'   Recommended values are \code{""}, \code{"-"}, or \code{"_"}.
#' @inheritParams texreg
#' @inheritParams matrixreg
#'
#' @author Philip Leifeld
#' @family texreg
#' @seealso \code{\link{texreg-package}} \code{\link{extract}}
#'
#' @references Leifeld, Philip (2013). texreg: Conversion of Statistical Model
#'   Output in R to LaTeX and HTML Tables. Journal of Statistical Software
#'   55(8): 1-24. \doi{10.18637/jss.v055.i08}.
#'
#' @examples
#' # Display models from ?lm:
#' ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
#' trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
#' group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
#' weight <- c(ctl, trt)
#' lm.D9 <- lm(weight ~ group)
#' lm.D90 <- lm(weight ~ group - 1)
#' screenreg(list(lm.D9, lm.D90))
#'
#' @export
screenreg <- function(l,
                      file = NULL,
                      single.row = FALSE,
                      stars = c(0.001, 0.01, 0.05),
                      custom.header = NULL,
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
  coltypes <- customcolumnnames(mod.names, custom.columns, custom.col.pos, types = TRUE)
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
  if (!is.character(outer.rule)) {
    stop("outer.rule must be a character.")
  } else if (nchar(outer.rule) > 1) {
    stop("outer.rule must be a character of maximum length 1.")
  } else if (outer.rule == "") {
    o.rule <- ""
  } else {
    o.rule <- paste(rep(outer.rule, table.width), collapse = "")
    string <- paste0(string, o.rule, "\n")
  }

  # specify multicolumn header
  spacing <- paste(rep(" ", column.spacing), collapse = "")
  mod.names <- c("", mod.names)
  if (!is.null(custom.header) && length(custom.header) > 0 && !any(is.na(custom.header))) {
    if (!"list" %in% class(custom.header) || length(custom.header) >= length(mod.names) || is.null(names(custom.header)) || !all(sapply(custom.header, is.numeric))) {
      stop("'custom.header' must be a named list of numeric vectors.")
    }
    ch <- unlist(custom.header)
    for (i in 1:length(ch)) {
      if (is.na(ch[i])) {
        stop("NA values are not permitted in 'custom.header'. Try leaving out the model indices that should not be included in the custom header.")
      }
      if (ch[i] %% 1 != 0) {
        stop("The model column indices in 'custom.header' must be provided as integer values.")
      }
      if (ch[i] < 1 || ch[i] >= length(mod.names)) {
        stop("The model column indices in 'custom.header' must be between 1 and the number of models.")
      }
      if (i > 1 && ch[i] <= ch[i - 1]) {
        stop("The model column indices in 'custom.header' must be strictly increasing.")
      }
    }

    ch <- rules <- paste0(rep(" ", nchar(output.matrix[1, 1])), collapse = "") # multicolumn header labels and mid-rules
    counter <- 0  # keeps track of how many columns we have processed to the left of the current start column
    for (i in 1:length(custom.header)) {
      # check if there are gaps within multicolumn blocks and throw error
      if (length(custom.header[[i]]) != custom.header[[i]][length(custom.header[[i]])] - custom.header[[i]][1] + 1) {
        stop("Each item in 'custom.header' must have strictly consecutive column indices, without gaps.")
      }

      # find out corrected column indices (ignoring the coef label column) of the current model after taking into account custom columns
      numCoefCol <- 0
      numCustomCol <- 0
      for (j in 1:length(coltypes)) {
        if (coltypes[j] == "coef") {
          numCoefCol <- numCoefCol + 1
        } else if (coltypes[j] == "customcol") {
          numCustomCol <- numCustomCol + 1
        }
        if (numCoefCol == custom.header[[i]][1]) {
          break() # break the loop if we have reached the number of models so far
        }
      }
      startIndex <- numCoefCol + numCustomCol # corrected column index
      numCoefCol <- 0
      numCustomCol <- 0
      for (j in 1:length(coltypes)) {
        if (coltypes[j] == "coef") {
          numCoefCol <- numCoefCol + 1
        } else if (coltypes[j] == "customcol") {
          numCustomCol <- numCustomCol + 1
        }
        if (numCoefCol == custom.header[[i]][length(custom.header[[i]])]) {
          break() # break the loop if we have reached the number of models so far
        }
      }
      stopIndex <- numCoefCol + numCustomCol # corrected column index

      # add empty cells for gaps and custom text columns
      if (startIndex > counter + 1) {
        spaces <- ""
        for (j in (startIndex):(counter + 2)) {
          spaces <- paste0(spaces, spacing, paste0(rep(" ", nchar(output.matrix[1, j])), collapse = ""))
          counter <- counter + 1
        }
        ch <- paste0(ch, spaces)
        rules <- paste0(rules, spaces)
      }

      # add multicolumn cells (+ 1 is for the coefficient label column)
      chars <- 0
      for (j in (startIndex + 1):(stopIndex + 1)) {
        chars <- chars + nchar(output.matrix[1, j]) + column.spacing
      }
      chars <- chars - column.spacing # last column does not have extra spacing
      label <- names(custom.header)[i]
      chars_label <- chars - nchar(label)
      if (chars_label < 0) {
        label <- substr(label, 1, chars)
      }
      before <- after <- (chars - nchar(label)) / 2
      if (before %% 1 != 0) {
        before <- before + 0.5
        after <- after - 0.5
      }
      ch <- paste0(ch,
                   spacing,
                   paste0(rep(" ", before), collapse = ""),
                   label,
                   paste0(rep(" ", after), collapse = "")
                   )
      rules <- paste0(rules,
                      spacing,
                      paste0(rep(inner.rule, chars), collapse = "")
                      )
      counter <- counter + stopIndex - startIndex + 1
    }
    string <- paste0(string, ch, "\n", rules, "\n")
  }

  # specify model names
  string <- paste(string, output.matrix[1, 1], sep = "")
  for (i in 2:ncol(output.matrix)) {
    string <- paste0(string, spacing, output.matrix[1, i])
  }
  string <- paste0(string, "\n")

  # mid rule 1
  if (!is.character(inner.rule)) {
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
  snote <- get_stars_note(stars = stars,
                          star.symbol = star.symbol,
                          symbol = symbol,
                          ci = ci,
                          ci.test = ci.test,
                          output = "ascii")

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

#' Convert regression output to a LaTeX table
#'
#' Conversion of \R regression output to a LaTeX table.
#'
#' The \code{\link{texreg}} function creates LaTeX code for inclusion in a LaTeX
#' document or for usage with \pkg{Sweave} or \pkg{knitr}, based on a list of
#' statistical models.
#'
#' @param custom.header An optional named list of multi-column headers that are
#'   placed above the model names. For example,
#'   \code{custom.header = list("abc" = 1:3, "ef" = 4:5)} will add the label
#'   \code{"abc"} to the first three models and \code{"ef"} to the fourth and
#'   fifth model. The column with coefficient names and any custom columns added
#'   by the \code{"custom.columns"} argument are not counted towards these
#'   positions. If \code{booktabs = TRUE}, \code{\\cmidrule} rules are added
#'   below the respective labels; otherwise \code{\\cline} lines are used.
#' @param file Using this argument, the resulting table is written to a file
#'   rather than to the \R prompt. The file name can be specified as a character
#'   string. Writing a table to a file can be useful for working with MS Office
#'   or LibreOffice. For example, using the \code{\link{htmlreg}} function, an
#'   HTML table can be written to a file with the extension \code{.doc} and
#'   opened with MS Word. The table can then be simply copied into any Word
#'   document, retaining the formatting of the table. Note that LibreOffice can
#'   import only plain HTML; CSS decorations are not supported; the resulting
#'   tables do not retain the full formatting in LibreOffice.
#' @param custom.note With this argument, a replacement text for the
#'   significance note below the table can be provided. If an empty
#'   \code{character} object is provided (\code{custom.note = ""}), the note
#'   will be omitted completely. If some character string is provided (e.g.,
#'   \code{custom.note = "My note"}), the significance legend is replaced by
#'   \code{My note}. The original significance legend can be included by
#'   inserting the \code{\%stars} wildcard. For example, a custom note can be
#'   added right after the significance legend by providing \code{custom.note =
#'   "\%stars. My note."}.
#'
#'   If the \code{threeparttable} argument is used, any note should be preceded
#'   by \code{"\\\\item"}, for example
#'   \code{"\\\\item \%stars. \\\\item Second note. \\\\item Third note."}, and
#'   it is possible to create line breaks in the formatted table by including
#'   \code{"\\\\\\\\"} and line breaks in the LaTeX code by including
#'   \code{"\\n"}, for example
#'   \code{"\\n\\\\item \%stars.\\\\\\\\\\n\\\\item Second line.\\n"}.
#' @param center Should the table be horizontally aligned at the center of the
#'   page?
#' @param caption Set the caption of the table.
#' @param caption.above Should the caption of the table be placed above the
#'   table? By default, it is placed below the table.
#' @param label Set the label of the \code{table} environment.
#' @param booktabs Use the \pkg{booktabs} LaTeX package to get thick horizontal
#'   rules in the output table (recommended).
#' @param lyx \code{logical}; if \code{TRUE}, each new line in the output is
#'   doubled, which facilitates transferring the output into the LyX document
#'   processor.
#' @param sideways If \code{sideways = TRUE} is set, the \code{table} floating
#'   environment is replaced by a \code{sidewaystable} float, and the
#'   \code{rotating} package is loaded in the preamble. The argument only has an
#'   effect if \code{table = TRUE} is also set.
#' @param longtable If \code{longtable = TRUE} is set, the \code{longtable}
#'   environment from the \code{longtable} LaTeX package is used to set tables
#'   across multiple pages. Note that this argument is not compatible with the
#'   \code{sideways} and \code{scalebox} arguments. These arguments will be
#'   automatically switched off when \code{longtable = TRUE} is set.
#' @param threeparttable If \code{threeparttable = TRUE} is set, the
#'   \code{threeparttable} environment will be used to enclose the
#'   \code{tabular} environment in the LaTeX code, and the significance note
#'   will be enclosed in a \code{tablenotes} environment. This permits word
#'   wrapping of long table notes and adequate spacing between multiple notes.
#'   See also the \code{custom.note} argument. If \code{longtable} is used,
#'   the \code{threeparttablex} LaTeX package is used instead of the
#'   \code{threeparttable} package.
#' @param use.packages If this argument is set to \code{TRUE} (= the default
#'   behavior), the required LaTeX packages are loaded in the beginning. If set
#'   to \code{FALSE}, the use package statements are omitted from the output.
#' @param table By default, \code{texreg} puts the actual \code{tabular} object
#'   in a \code{table} floating environment. To get only the \code{tabular}
#'   object without the whole table header, set \code{table = FALSE}.
#' @param tabular By default, the table contents are wrapped in a \code{tabular}
#'   environment. To get only the contents for each row without the environment,
#'   set \code{tabular = FALSE}. Note that if \code{tabular = FALSE}, the
#'   \code{table} argument must also be \code{FALSE}, otherwise a warning is
#'   printed. Switching off the tabular environment may be useful for designing
#'   one's own table more flexibly, for example using \code{tabular*} or
#'   \code{tabularx} environments in LaTeX.
#' @param no.margin In order to save space, inner margins of tables can be
#'   switched off.
#' @param fontsize The \code{fontsize} argument serves to change the font size
#'   used in the table. Valid values are \code{"tiny"}, \code{"scriptsize"},
#'   \code{"footnotesize"}, \code{"small"}, \code{"normalsize"}, \code{"large"},
#'   \code{"Large"}, \code{"LARGE"}, \code{"huge"}, and \code{"Huge"}. Note that
#'   the \code{scalebox} argument often achieves better results when the goal is
#'   to change the size of the table.
#' @param scalebox The \code{scalebox} argument serves to resize the table. For
#'   example, \code{scalebox = 1.0} is equivalent to the normal size,
#'   \code{scalebox = 0.5} decreases the size of the table by one half, and
#'   \code{scalebox = 2.0} doubles the space occupied by the table. Note that
#'   the \code{scalebox} argument does not work when the \code{longtable}
#'   argument is used.
#' @param float.pos This argument specifies where the table should be located on
#'   the page or in the document. By default, no floating position is specified,
#'   and LaTeX takes care of the position automatically. Possible values include
#'   \code{"h"} (here), \code{"p"} (page), \code{"t"} (top), \code{"b"}
#'   (bottom), any combination thereof, e.g., \code{"tb"}, or any of these
#'   values followed by an exclamation mark, e.g. \code{"t!"}, in order to
#'   enforce this position. The square brackets do not have to be specified.
#' @return A \code{character} object with a regression table and LaTeX markup.
#'   The object has an additional \code{"texregTable"} class identifier, which
#'   causes the object to be formatted nicely on screen when printed.
#' @inheritParams matrixreg
#'
#' @author Philip Leifeld
#' @family texreg
#' @seealso \code{\link{texreg-package}} \code{\link{extract}}
#'
#' @references Leifeld, Philip (2013). texreg: Conversion of Statistical Model
#'   Output in R to LaTeX and HTML Tables. Journal of Statistical Software
#'   55(8): 1-24. \doi{10.18637/jss.v055.i08}.
#'
#' @keywords print misc utilities IO programming|interface
#'
#' @examples
#' # Linear mixed-effects models
#' library("nlme")
#' model.1 <- lme(distance ~ age, data = Orthodont, random = ~ 1)
#' model.2 <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1)
#' texreg(list(model.1, model.2), booktabs = TRUE, dcolumn = TRUE)
#'
#' # Ordinary least squares model (example from the 'lm' help file)
#' ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#' trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#' group <- gl(2,10,20, labels = c("Ctl","Trt"))
#' weight <- c(ctl, trt)
#' lm.D9 <- lm(weight ~ group)
#' table.string <- texreg(lm.D9, return.string = TRUE)
#' cat(table.string)
#'
#' @export
texreg <- function(l,
                   file = NULL,
                   single.row = FALSE,
                   stars = c(0.001, 0.01, 0.05),
                   custom.header = NULL,
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
                   siunitx = FALSE,
                   lyx = FALSE,
                   sideways = FALSE,
                   longtable = FALSE,
                   threeparttable = FALSE,
                   use.packages = TRUE,
                   table = TRUE,
                   tabular = TRUE,
                   no.margin = FALSE,
                   fontsize = NULL,
                   scalebox = NULL,
                   float.pos = "",
                   ...) {

  # check dcolumn vs. siunitx
  if (isTRUE(dcolumn) && isTRUE(siunitx)) {
    dcolumn <- FALSE
    msg <- paste("The dcolumn and siunitx packages cannot be used at",
                 "the same time. Switching off 'dcolumn'.")
    warning(msg)
  }

  # check dcolumn vs. bold
  if (isTRUE(dcolumn) && bold > 0) {
    dcolumn <- FALSE
    msg <- paste("The dcolumn package and the 'bold' argument cannot be used at",
                 "the same time. Switching off 'dcolumn'.")
    if (length(stars) > 1 || (length(stars) > 0 && (stars == TRUE || stars != 0))) {
      warning(paste(msg, "You should also consider setting stars = 0."))
    } else {
      warning(msg)
    }
  }

  # check siunitx vs. bold
  if (isTRUE(siunitx) && bold > 0) {
    siunitx <- FALSE
    msg <- paste("The siunitx package and the 'bold' argument cannot be used at",
                 "the same time. Switching off 'siunitx'.")
    if (length(stars) > 1 || (length(stars) > 0 && (stars == TRUE || stars != 0))) {
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

  # check longtable vs. scalebox
  if (isTRUE(longtable) && !is.null(scalebox)) {
    scalebox <- NULL
    warning(paste("'longtable' and 'scalebox' are not compatible. Setting scalebox = NULL."))
  }

  # check table vs. tabular
  if (isTRUE(table) && !isTRUE(tabular)) {
    table <- FALSE
    warning("Setting 'table = FALSE' because 'tabular = FALSE'.")
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
                             siunitx = siunitx,
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
    if (isTRUE(single.row) && coltypes[i] == "coef" && !isTRUE(siunitx)) {
      if (isTRUE(ci[coefcount])) {
        separator <- "]"
      } else {
        separator <- ")"
      }
    } else {
      separator <- "."
    }
    if (coltypes[i] %in% c("coef", "customcol")) {
      alignmentletter <- "c"
    } else if (coltypes[i] == "coefnames") {
      alignmentletter <- "l"
    }
    if (isTRUE(dcolumn)) {
      if (coltypes[i] != "coef") {
        coldef <- paste0(coldef, alignmentletter, margin.arg, " ")
      } else {
        dl <- compute.width(output.matrix[, i],
                            left = TRUE,
                            single.row = single.row,
                            bracket = separator)
        dr <- compute.width(output.matrix[, i],
                            left = FALSE,
                            single.row = single.row,
                            bracket = separator)
        coldef <- paste0(coldef, "D{", separator, "}{", separator, "}{",
                         dl, separator, dr, "}", margin.arg, " ")
      }
    } else if (isTRUE(siunitx)) {
      if (coltypes[i] != "coef") {
        coldef <- paste0(coldef, alignmentletter, margin.arg, " ")
      } else {
        dl <- compute.width(output.matrix[, i],
                            left = TRUE,
                            single.row = single.row,
                            bracket = separator)
        dr <- compute.width(output.matrix[, i],
                            left = FALSE,
                            single.row = single.row,
                            bracket = separator)
        coldef <- paste0(coldef, "S[table-format=", dl, separator, dr, "]", margin.arg, " ")
      }
    } else {
      coldef <- paste0(coldef, alignmentletter, margin.arg, " ")
    }
  }
  coldef <- trimws(coldef) # remove white space at the end due to the for-loop

  string <- "\n"
  linesep <- ifelse(isTRUE(lyx), "\n\n", "\n")

  # write table header
  if (isTRUE(use.packages)) {
    if (!is.null(scalebox)) {
      string <- paste0(string, "\\usepackage{graphicx}", linesep)
    }
    if (isTRUE(sideways) & isTRUE(table)) {
      string <- paste0(string, "\\usepackage{rotating}", linesep)
    }
    if (isTRUE(booktabs)) {
      string <- paste0(string, "\\usepackage{booktabs}", linesep)
    }
    if (isTRUE(dcolumn)) {
      string <- paste0(string, "\\usepackage{dcolumn}", linesep)
    }
    if (isTRUE(siunitx)) {
      string <- paste0(string, "\\usepackage{siunitx}", linesep)
    }
    if (isTRUE(longtable)) {
      string <- paste0(string, "\\usepackage{longtable}", linesep)
    }
    if (isTRUE(threeparttable)) {
      if (isTRUE(longtable)) {
        string <- paste0(string, "\\usepackage{threeparttablex}", linesep)
      } else {
        string <- paste0(string, "\\usepackage{threeparttable}", linesep)
      }
    }
    if (!is.null(scalebox) || isTRUE(dcolumn) || isTRUE(siunitx) || isTRUE(booktabs) || isTRUE(sideways) || isTRUE(longtable) || isTRUE(threeparttable)) {
      string <- paste0(string, linesep)
    }
  }

  # stars note (define now, add later)
  snote <- get_stars_note(stars = stars,
                          star.symbol = "*",
                          symbol = symbol,
                          ci = ci,
                          ci.test = ci.test,
                          output = "latex")

  if (is.null(fontsize)) {
    notesize <- "scriptsize"
  } else if (fontsize == "tiny" ||
             fontsize == "scriptsize" ||
             fontsize == "footnotesize" ||
             fontsize == "small") {
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
    } else if (!isTRUE(threeparttable)) {
      note <- paste0("\\multicolumn{", length(mod.names),
                     "}{l}{\\", notesize, "{", snote, "}}")
    } else {
      note <- paste0("\\", notesize, "{\\item ", snote, "}")
    }
  } else if (custom.note == "") {
    note <- ""
  } else {
    if (!isTRUE(threeparttable)) {
      note <- paste0("\\multicolumn{", length(mod.names),
                     "}{l}{\\", notesize, "{", custom.note, "}}")
    } else {
      note <- paste0("\\", notesize, "{", custom.note, "}")
    }
    note <- gsub("%stars", snote, note, perl = TRUE)
  }
  if (note != "") {
    if (isTRUE(longtable) && !isTRUE(threeparttable)) {  # longtable requires line break after note/caption
      note <- paste0(note, "\\\\", linesep)
    } else {
      note <- paste0(note, linesep)
    }
  }

  if (isTRUE(longtable)) {
    if (isTRUE(center)) {
      string <- paste0(string, "\\begin{center}\n")
    }
    if (!is.null(fontsize)) {
      string <- paste0(string, "\\begin{", fontsize, "}", linesep)
    }
    if (isTRUE(threeparttable)) {
      string <- paste0(string,
                       "\\begin{ThreePartTable}",
                       linesep,
                       "\\begin{TableNotes}[flushleft]",
                       linesep,
                       note,
                       "\\end{TableNotes}",
                       linesep)
    }
    if (float.pos == "") {
      string <- paste0(string, "\\begin{longtable}{", coldef, "}", linesep)
    } else {
      string <- paste0(string, "\\begin{longtable}[", float.pos, "]{", coldef, "}", linesep)
    }
  } else {  # table or sidewaystable
    if (isTRUE(table)) {
      if (isTRUE(sideways)) {
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
      if (isTRUE(caption.above)) {
        string <- paste0(string, "\\caption{", caption, "}", linesep)
      }
      if (isTRUE(center)) {
        string <- paste0(string, "\\begin{center}", linesep)
      }
      if (!is.null(fontsize)) {
        string <- paste0(string, "\\begin{", fontsize, "}", linesep)
      }
      if (!is.null(scalebox)) {
        string <- paste0(string, "\\scalebox{", scalebox, "}{", linesep)
      }
    }
    if (isTRUE(siunitx)) {
      string <- paste0(string, "\\sisetup{parse-numbers=false, table-text-alignment=right}", linesep)
    }
    if (isTRUE(threeparttable)) {
      string <- paste0(string, "\\begin{threeparttable}", linesep)
    }
    if (isTRUE(tabular)) {
      string <- paste0(string, "\\begin{tabular}{", coldef, "}", linesep)
    }
  }

  # horizontal rule above the table
  tablehead <- ""
  if (isTRUE(booktabs)) {
    tablehead <- paste0(tablehead, "\\toprule", linesep)
  } else {
    tablehead <- paste0(tablehead, "\\hline", linesep)
  }

  # specify multicolumn header
  if (!is.null(custom.header) && length(custom.header) > 0 && !any(is.na(custom.header))) {
    if (!"list" %in% class(custom.header) || length(custom.header) >= length(mod.names) || is.null(names(custom.header)) || !all(sapply(custom.header, is.numeric))) {
      stop("'custom.header' must be a named list of numeric vectors.")
    }
    ch <- unlist(custom.header)
    for (i in 1:length(ch)) {
      if (is.na(ch[i])) {
        stop("NA values are not permitted in 'custom.header'. Try leaving out the model indices that should not be included in the custom header.")
      }
      if (ch[i] %% 1 != 0) {
        stop("The model column indices in 'custom.header' must be provided as integer values.")
      }
      if (ch[i] < 1 || ch[i] >= length(mod.names)) {
          stop("The model column indices in 'custom.header' must be between 1 and the number of models.")
      }
      if (i > 1 && ch[i] <= ch[i - 1]) {
          stop("The model column indices in 'custom.header' must be strictly increasing.")
      }
    }
    ch <- ""      # multicolumn header labels
    rules <- ""   # multicolumn header mid-rules (booktabs package if possible)
    counter <- 0  # keeps track of how many columns we have processed to the left of the current start column
    for (i in 1:length(custom.header)) {
      # check if there are gaps within multicolumn blocks and throw error
      if (length(custom.header[[i]]) != custom.header[[i]][length(custom.header[[i]])] - custom.header[[i]][1] + 1) {
        stop("Each item in 'custom.header' must have strictly consecutive column indices, without gaps.")
      }

      # find out corrected column indices (ignoring the coef label column) of the current model after taking into account custom columns
      numCoefCol <- 0
      numCustomCol <- 0
      for (j in 1:length(coltypes)) {
        if (coltypes[j] == "coef") {
          numCoefCol <- numCoefCol + 1
        } else if (coltypes[j] == "customcol") {
          numCustomCol <- numCustomCol + 1
        }
        if (numCoefCol == custom.header[[i]][1]) {
          break() # break the loop if we have reached the number of models so far
        }
      }
      startIndex <- numCoefCol + numCustomCol # corrected column index
      numCoefCol <- 0
      numCustomCol <- 0
      for (j in 1:length(coltypes)) {
        if (coltypes[j] == "coef") {
          numCoefCol <- numCoefCol + 1
        } else if (coltypes[j] == "customcol") {
          numCustomCol <- numCustomCol + 1
        }
        if (numCoefCol == custom.header[[i]][length(custom.header[[i]])]) {
          break() # break the loop if we have reached the number of models so far
        }
      }
      stopIndex <- numCoefCol + numCustomCol # corrected column index

      # check if there are gaps between multicolumn headers and fill up with empty cells
      if (i > 1) {
        emptycells <- custom.header[[i]][1] - custom.header[[i - 1]][length(custom.header[[i - 1]])] - 1
        if (emptycells > 0) {
          for (j in 1:emptycells) {
            ch <- paste0(ch, " &")
            counter <- counter + 1
          }
        }
      }

      # add empty cells also for custom text columns
      if (startIndex > counter + 1) {
        difference <- startIndex - (counter + 1)
        for (j in 1:difference) {
          ch <- paste0(ch, " &")
          counter <- counter + 1
        }
      }

      # add multicolumn cells (+ 1 is for the coefficient label column)
      ch <- paste0(ch, " & \\multicolumn{", stopIndex - startIndex + 1, "}{c}{", names(custom.header)[i], "}")
      counter <- counter + stopIndex - startIndex + 1

      # add mid-rules (+ 1 is for the coefficient label column)
      if (isTRUE(booktabs)) {
        rules <- paste0(rules, ifelse(i > 1, " ", ""), "\\cmidrule(lr){", startIndex + 1, "-", stopIndex + 1, "}")
      } else {
        rules <- paste0(rules, ifelse(i > 1, " ", ""), "\\cline{", startIndex + 1, "-", stopIndex + 1, "}")
      }
    }
    tablehead <- paste0(tablehead, ch, " \\\\", linesep, rules, linesep)
  }

  # specify model names
  tablehead <- paste0(tablehead, mod.names[1])
  if (isTRUE(dcolumn)) {
    for (i in 2:length(mod.names)) {
      if (coltypes[i] != "coef") {
        tablehead <- paste0(tablehead, " & ", mod.names[i])
      } else {
        tablehead <- paste0(tablehead, " & \\multicolumn{1}{c}{", mod.names[i], "}")
      }
    }
  } else if (isTRUE(siunitx)) {
    for (i in 2:length(mod.names)) {
      if (coltypes[i] != "coef") {
        tablehead <- paste0(tablehead, " & ", mod.names[i])
      } else {
        tablehead <- paste0(tablehead, " & {", mod.names[i], "}")
      }
    }
  } else {
    for (i in 2:length(mod.names)) {
      tablehead <- paste0(tablehead, " & ", mod.names[i])
    }
  }

  # horizontal rule between model names and coefficients (define now, add later)
  if (isTRUE(booktabs)) {
    tablehead <- paste0(tablehead, " \\\\", linesep, "\\midrule", linesep)
  } else {
    tablehead <- paste0(tablehead, " \\\\", linesep, "\\hline", linesep)
  }
  if (isFALSE(longtable)) {
    string <- paste0(string, tablehead)
  }

  # bottom rule (define now, add later)
  if (isTRUE(booktabs)) {
    bottomline <- paste0("\\bottomrule", linesep)
  } else {
    bottomline <- paste0("\\hline", linesep)
  }

  # write table header (and footer, in the case of longtable)
  if (isTRUE(longtable)) {
    if (isTRUE(caption.above)) {
      string <- paste0(string, "\\caption{", caption, "}", linesep, "\\label{",
                       label, "}\\\\", linesep, tablehead, "\\endfirsthead", linesep,
                       tablehead, "\\endhead", linesep, bottomline, "\\endfoot", linesep,
                       bottomline, ifelse(isTRUE(threeparttable), "\\insertTableNotes\\\\\n", note),
                       "\\endlastfoot", linesep)
    } else {
      string <- paste0(string, tablehead, "\\endfirsthead", linesep, tablehead,
                       "\\endhead", linesep, bottomline, "\\endfoot", linesep, bottomline,
                       ifelse(isTRUE(threeparttable), "\\insertTableNotes\\\\\n", note),
                       "\\caption{", caption, "}", linesep, "\\label{", label, "}",
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
    if (isTRUE(booktabs)) {
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
  if (isFALSE(longtable)) {
    string <- paste0(string, bottomline)
    if (isTRUE(threeparttable)) {
      string <- paste0(string,
                       ifelse(isTRUE(tabular), "\\end{tabular}", ""),
                       ifelse(isTRUE(tabular), linesep, ""),
                       "\\begin{tablenotes}[flushleft]",
                       linesep,
                       note,
                       "\\end{tablenotes}",
                       linesep,
                       "\\end{threeparttable}",
                       linesep
                       )
    } else {
      if(isTRUE(tabular)) {
        string <- paste0(string, note, "\\end{tabular}", linesep)
      } else {
        string <- paste0(string, note)
      }
    }
  }

  # take care of center, scalebox and table environment
  if (isTRUE(longtable)) {
    string <- paste0(string, "\\end{longtable}", linesep)
    if (isTRUE(threeparttable)) {
      string <- paste0(string, "\\end{ThreePartTable}", linesep)
    }
    if (!is.null(fontsize)) {
      string <- paste0(string, "\\end{", fontsize, "}", linesep)
    }
    if (isTRUE(center)) {
      string <- paste0(string, "\\end{center}", linesep)
    }
  } else if (isTRUE(table)) {
    if (!is.null(fontsize)) {
      string <- paste0(string, "\\end{", fontsize, "}", linesep)
    }
    if (!is.null(scalebox)) {
      string <- paste0(string, "}", linesep)
    }
    if (isFALSE(caption.above)) {
      string <- paste0(string, "\\caption{", caption, "}", linesep)
    }
    string <- paste0(string, "\\label{", label, "}", linesep)
    if (isTRUE(center)) {
      string <- paste0(string, "\\end{center}", linesep)
    }
    if (isTRUE(sideways)) {
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

#' Export regression output to an MS Word file
#'
#' Export regression output to an MS Word file.
#'
#' The \code{wordreg} function creates a Microsoft Word document with the
#' requested table.
#'
#' @inheritParams matrixreg
#' @inheritParams texreg
#'
#' @author Vincent Arel-Bundock
#' @family texreg
#' @seealso \code{\link{texreg-package}} \code{\link{extract}}
#'
#' @examples
#' \dontrun{
#' # Use models from ?lm:
#' ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
#' trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
#' group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
#' weight <- c(ctl, trt)
#' lm.D9 <- lm(weight ~ group)
#' lm.D90 <- lm(weight ~ group - 1)
#' wordreg(list(lm.D9, lm.D90), file = "testfile.doc")
#' unlink("testfile.doc")
#' }
#'
#' @export
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
                    ...) {

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
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
                   output.type = "ascii",
                   include.attributes = FALSE,
                   trim = TRUE,
                   ...
  )
  wd <- getwd()
  f = tempfile(fileext = ".Rmd")
  if (!all(mat[1, ] == "")) { # column names for kable
    dvnames <- mat[1, ]
    mat <- mat[-1, ]
  } else {
    dvnames <- NA
  }
  cat(file = f, "```{r, echo = FALSE}
                    knitr::kable(mat, col.names = dvnames)
                    ```", append = TRUE)
  rmarkdown::render(f, output_file = paste0(wd, "/", file))
}


# Internal helpers -------------------------------------------------------------

#' Display version number and date when the package is loaded.
#' @importFrom utils packageDescription
#' @noRd
.onAttach <- function(libname, pkgname) {
  desc  <- utils::packageDescription(pkgname, libname)
  packageStartupMessage(
    "Version:  ", desc$Version, "\n",
    "Date:     ", desc$Date, "\n",
    "Author:   ", "Philip Leifeld (University of Manchester)", "\n\n",
    "Consider submitting praise using the praise or praise_interactive functions.\n",
    "Please cite the JSS article in your publications -- see citation(\"texreg\")."
  )
}

#' Convert a number into a string with rounded decimal places
#'
#' Reformat a coefficient as a string with a certain number of decimal places.
#'
#' This function takes a \code{numeric} object, usually a coefficient from a
#' statistical model, and converts it into a \code{character} object. The user
#' can choose to how many decimal places the number is rounded (usually two in
#' most published regression models) and whether there should be a leading zero
#' if the coefficient is between 0 and 1.
#'
#' @param x A \code{numeric} object to reformat.
#' @param lead.zero Should a leading zero be printed if the coefficient is
#'   non-negative and smaller than one (e.g., \code{"0.23"} as opposed to
#'   \code{".23"})?
#' @param digits Number of decimal places to round to.
#' @return A reformatted coefficient string as a \code{character} object.
#'
#' @keywords internal
#' @author Philip Leifeld
#' @seealso \code{\link{texreg}}
#'
#' @export
coeftostring <- function(x, lead.zero = FALSE, digits = 2) {
  if (is.character(x)) {
    return(x)
  }
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

#' Compute maximum column width left and right of a decimal separator
#'
#' Compute maximum column width left and right of a decimal separator.
#'
#' This function takes a vector of \code{character} objects with coefficients,
#' usually a column of a regression table, and computes the maximal width left
#' or right of the decimal separator or bracket at which the cells are aligned
#' vertically. This is useful in the context of the \code{\link{texreg}}
#' function when the \code{dcolumn} or \code{siunitx} arguments are used for
#' vertical decimal point alignment.
#'
#' @param v A \code{character} vector representing a column in a regression
#'   table.
#' @param left Should the width left of the separator/bracket be calculated? If
#'   \code{FALSE}, the width right of the separator/bracket is computed.
#' @param single.row Was the \code{single.row} argument used to construct the
#'   regression table? I.e., are both the coefficient and uncertainty measure
#'   (SE or CI) in the same rows of the matrix?
#' @param bracket The separator symbol to match. These can be closing
#'   parentheses (in the case of standard errors when \code{single.row} is
#'   switched on), closing square brackets (in the case of confidence
#'   intervals), or dots (in the case of \code{single.row = FALSE}, for decimal
#'   alignment at the actual decimal separator).
#' @return A number indicating the maximal width left or right of the separator.
#'
#' @keywords internal
#' @author Philip Leifeld
#' @seealso \code{\link{texreg}}
#'
#' @export
compute.width <- function(v, left = TRUE, single.row = FALSE, bracket = ")") {
  if (isFALSE(single.row)) {
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
        if (length(ssp[[i]]) > 1) {
          r <- ""
          for (j in 2:length(ssp[[i]])) {
            r <- paste(r, ssp[[i]][j], sep = ".")
          }
          right.side <- append(right.side, r)
        }
      }
    }
  }
  if (isTRUE(left)) {
    left.side <- sub(" \\\\; ", " ", left.side)
    left.side <- gsub("[.]", "", left.side) # do not count decimal separators on RHS because they are quite narrow
    # correct for custom.gof.rows with text, which is set as \multicolumn or {}
    left.side <- left.side[!grepl("\\\\multicolumn[{]1[}][{]c[}][{]", left.side)]
    left.side <- left.side[!grepl("^\\{", left.side)]
    v.length <- max(c(0, nchar(left.side)), na.rm = TRUE)
  } else {
    right.side <- sub("\\^\\{", "", right.side) # only count stars, not the ^{} around them
    right.side <- sub("\\}", "", right.side)
    right.side <- gsub(" \\\\; ", " ", right.side) # replace guaranteed space by normal space for counting purposes
    right.side <- gsub("[.]", "", right.side) # do not count decimal separators on RHS because they are quite narrow
    v.length <- max(c(0, nchar(right.side)), na.rm = TRUE)
  }
  return(v.length)
}

#' Determine column names or column types if custom columns are present
#'
#' Determine column names or column types if custom columns are present.
#'
#' This function takes model names (as saved in the attributes of a matrix
#' generated by \code{\link{matrixreg}}, for example) and the
#' \code{custom.columns} and \code{custom.col.pos} arguments of
#' \code{link{matrixreg}} or related functions and determines the column types
#' (\code{"coefnames"}, \code{"coef"}, or \code{"customcol"}) or model names in
#' the presence of custom columns.
#'
#' @param modelnames A \code{character} vector of original model names, before
#'   any custom columns are inserted.
#' @param custom.columns The same argument as specified in the
#'   \code{link{matrixreg}} function.
#' @param custom.col.pos The same argument as specified in the
#'   \code{link{matrixreg}} function.
#' @param types Return the column types? If \code{FALSE}, the column names in
#'   the possible presence of custom columns are returned.
#' @return A \code{character} vector with column names or types in the possible
#'   presence of custom columns. If \code{types = TRUE}, the vector contains
#'   the values \code{"coefnames"} (for the first column), \code{"coef"} (for
#'   columns with coefficients), or \code{"customcol"} (for custom new columns).
#'
#' @keywords internal
#' @author Philip Leifeld
#' @seealso \code{\link{matrixreg}}
#'
#' @export
customcolumnnames <- function(modelnames,
                              custom.columns,
                              custom.col.pos,
                              types = FALSE) {

  # adjust arguments
  modelnames <- c("", modelnames)
  if (is.null(custom.columns)) {
    if (isFALSE(types)) {
      return(modelnames)
    } else {
      return(c("coefnames", rep("coef", length(modelnames) - 1)))
    }
  }
  if (!"list" %in% class(custom.columns)) {
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

#' Extract all data necessary for generating a table from statistical models
#'
#' Extract all data necessary for generating a table from a list of models.
#'
#' This function applies the \code{link{extract}} function and its respective
#' methods to each element in a list of statistical models in order to extract
#' coefficients, standard errors, p-values, confidence intervals, and
#' goodness-of-fit statistics for generating a regression table.
#'
#' @param l A list of statistical models.
#' @param ... Arguments to be passed over to the \code{\link{extract}} function,
#'   such as \code{include.nobs} or \code{include.aic}. Details are provided in
#'   the documentation of the \code{\link{extract}} function.
#' @return A list of \linkS4class{texreg} objects.
#'
#' @keywords internal
#' @author Philip Leifeld
#' @seealso \code{\link{extract}}
#'
#' @export
get.data <- function(l, ...) {
  # if a single model is handed over, put model inside a list
  if (!"list" %in% class(l)[1]) {
    l <- list(l)
  }

  # create list of texreg objects
  models <- NULL
  for (i in 1:length(l)) {
    model <- extract(l[[i]], ...)
    if ("list" %in% class(model)) { # must be a nested list of models (e.g., systemfit)
      models <- append(models, model)
    } else { # normal case; one model
      models <- append(models, list(model))
    }
  }
  return(models)
}

#' Create a legend for the stars in a regression table
#'
#' Create a legend for the stars in a regression table.
#'
#' This function creates a stars note as a legend to be placed below a
#' regression table. The note contains the p-value or confidence interval
#' significance levels and stars attached to them.
#'
#' @param stars A numeric vector of cut-offs, with a maximum of four numbers.
#' @param star.symbol The character to repeat for the first three levels of
#'   significance.
#' @param symbol The character for the fourth level of significance.
#' @param ci Confidence intervals instead of standard errors?
#' @param ci.test The null hypothesis value, for example \code{0} (the normal
#'   case) or \code{1} (e.g., with exponentiated coefficients). A star is added
#'   if this value is outside the confidence interval.
#' @param output The output type of the note. This can be \code{"ascii"},
#'   \code{"latex"}, or \code{"html"}.
#' @return A \code{character} string to be put below the regression table. It
#'   describes the thresholds for the significance stars.
#'
#' @keywords internal
#' @author Philip Leifeld
#' @seealso \code{\link{texreg}}, \code{\link{htmlreg}}, \code{\link{screenreg}}
#'
#' @export
get_stars_note <- function(stars = c(0.01, 0.05, 0.1),
                           star.symbol = "*",
                           symbol = ".",
                           ci = FALSE,
                           ci.test = NULL,
                           output = "ascii") {

  # sanity checks and prep
  if (!output %in% c("ascii", "latex", "html")) {
    stop("'output' argument must be 'ascii', 'latex', or 'html'.")
  }
  if (!is.numeric(ci.test) && !is.null(ci.test) && !any(is.na(ci.test))) {
    stop("The argument 'ci.test' must be NULL, NA, or numeric.")
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
      p_note <- paste0(symbols,
                       " p < ", st)
    } else if (output == "latex") {
      p_note <- paste0("$^{",
                       symbols,
                       "}p<",
                       st,
                       "$")
    } else if (output == "html") {
      p_note <- paste0("<sup>",
                       symbols,
                       "</sup>p &lt; ",
                       st)
    }
    p_note <- paste(p_note, collapse = "; ")
  } else { # no stars supplied = empty note
    p_note <- ""
  }

  # ci_note
  if (ci_note_flag) {  # ci calculated for at least one model -> build ci note
    if (is.numeric(ci.test) && all(!is.na(ci.test))) { # sanity check
      if (length(ci.test) == 1) {
        ci_note <- paste(ci.test, "outside the confidence interval")
      } else {
        ci_note <- "Null hypothesis value outside the confidence interval"
      }
    } else {
      ci_note <- ""
    }
    if (output == "ascii") {
      ci_symbol <- "*"
    } else if (output == "latex") {
      ci_symbol <- "$^*$"
    } else if (output == "html") {
      ci_symbol <- paste0("<sup>&#42;</sup>")
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

  return(s_note)
}

#' Replace symbols in a character string or vector by LaTeX equivalents
#'
#' Replace symbols in a character string or vector by LaTeX equivalents-
#'
#' This function is an internal helper function that takes a \code{character}
#' object or vector and replaces symbols like underscores, angle brackets, or
#' superscripted numbers by properly escaped LaTeX equivalents in order not to
#' break any LaTeX table code when the \code{link{texreg}} function is used.
#'
#' @param x A \code{character} object of arbitrary length.
#' @return Same as the input object but with escaped or replaced symbols.
#'
#' @keywords internal
#' @author Philip Leifeld
#' @seealso \code{\link{texreg}}
#'
#' @export
names2latex <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  x <- sapply(x, function(a) {
    if (!is.na(a) && !grepl("\\$", a)) {
      a <- gsub("_", "\\\\_", a)
      a <- gsub("<", "\\$<\\$", a)
      a <- gsub(">", "\\$>\\$", a)
      a <- gsub("%", "\\\\%", a)
      a <- gsub("\\^2", "\\$^2\\$", a)
      a <- gsub("\\^3", "\\$^3\\$", a)
      a <- gsub("\\^4", "\\$^4\\$", a)
      a <- gsub("\\^5", "\\$^5\\$", a)
    }
    return(a)
  })
  return(x)
}

#' Replace coefs, SEs, p-values, and/or CIs by custom values if provided
#'
#' Replace coefs, SEs, p-values, and/or CIs by custom values if provided.
#'
#' This function replaces coefficients, standard errors, p-values, and/or
#' confidence intervals in a list of \linkS4class{texreg} objects. It is used by
#' the \code{\link{matrixreg}} and \code{\link{plotreg}} functions. The new
#' values must be provided as lists of equal length as the list of models, with
#' each element representing a vector of replacement values. If the arguments
#' have value \code{0}, the original values are retained. More details are found
#' in the documentation of the \code{\link{matrixreg}} function.
#'
#' @param models A list of \linkS4class{texreg} objects (e.g., as returned by
#'   \code{\link{get.data}}).
#' @param override.coef Replacement list of coefficient vectors.
#' @param override.se Replacement list of standard error vectors
#' @param override.pvalues Replacement list of p-value vectors.
#' @param override.ci.low Replacement list of lower-bound confidence interval
#'   values.
#' @param override.ci.up Replacement list of upper-bound confidence interval
#'   values.
#' @return Same list as input list of models, but with replaced values.
#'
#' @keywords internal
#' @author Philip Leifeld
#' @seealso \code{\link{matrixreg}}, \code{\link{plotreg}}, \code{\link{texreg}}
#'
#' @export
override <- function(models,
                     override.coef = 0,
                     override.se = 0,
                     override.pvalues = 0,
                     override.ci.low = 0,
                     override.ci.up = 0) {
  # check validity of override arguments for p-values and SEs
  if ("list" %in% class(override.se) || length(override.se) > 1 || override.se[1] != 0) {
    if (length(override.pvalues) == 1 && !"list" %in% class(override.pvalues) && override.pvalues[1] == 0) {
      warning("Standard errors were provided using 'override.se', but p-values were not replaced!")
    }
  }
  if ("list" %in% class(override.pvalues) || length(override.pvalues) > 1 || override.pvalues[1] != 0) {
    if (length(override.se) == 1 && !"list" %in% class(override.se) && override.se[1] == 0) {
      warning("p-values were provided using 'override.pvalues', but standard errors were not replaced!")
    }
  }

  # replace coefs, SEs, p-values, and/or CIs by custom values if provided
  for (i in 1:length(models)) {
    # override coefficients
    if (!"list" %in% class(override.coef) && length(override.coef) == 1 && override.coef == 0) {
      cf <- models[[i]]@coef
    } else if (is.numeric(override.coef) &&
               length(models) == 1 &&
               length(override.coef) == length(models[[i]]@coef)) {
      cf <- override.coef
    } else if (!"list" %in% class(override.coef)) {
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
    } else if (!is.numeric(override.coef[[i]])) {
      warning("Coefficients provided for model", i, "are not numeric. Using default values.")
      cf <- models[[i]]@coef
    } else {
      cf <- override.coef[[i]]
    }
    models[[i]]@coef <- cf

    # override standard errors
    if (!"list" %in% class(override.se) && length(override.se) == 1 && override.se == 0) {
      se <- models[[i]]@se
    } else if (is.numeric(override.se) &&
               length(models) == 1 &&
               length(override.se) == length(models[[i]]@se)) {
      se <- override.se
    } else if (!"list" %in% class(override.se)) {
      warning("SEs must be provided as a list. Using default SEs.")
      se <- models[[i]]@se
    } else if (length(override.se) != length(models)) {
      warning("Number of SEs provided does not match number of models. Using default SEs.")
      se <- models[[i]]@se
    } else if (length(models[[i]]@se) != length(override.se[[i]])) {
      warning(paste0("Number of SEs provided does not match number of ",
                     "coefficients in model ", i, ". Using default SEs."))
      se <- models[[i]]@se
    } else if (!is.numeric(override.se[[i]])) {
      warning(paste("SEs provided for model", i, "are not numeric. Using default SEs."))
      se <- models[[i]]@se
    } else {
      se <- override.se[[i]]
    }
    models[[i]]@se <- se

    # override p-values
    if (!"list" %in% class(override.pvalues) && length(override.pvalues) == 1 && override.pvalues == 0) {
      pval <- models[[i]]@pvalues
    } else if (is.numeric(override.pvalues) &&
               length(models) == 1 &&
               length(override.pvalues) == length(models[[i]]@pvalues)) {
      pval <- override.pvalues
    } else if (!"list" %in% class(override.pvalues)) {
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
    } else if (!is.numeric(override.pvalues[[i]])) {
      warning(paste("p-values provided for model", i, "are not numeric. Using default p-values."))
      pval <- models[[i]]@pvalues
    } else {
      pval <- override.pvalues[[i]]
    }
    models[[i]]@pvalues <- pval

    # override lower bound of confidence intervals
    if (is.null(override.ci.low)) {
      # do nothing
    } else if (!"list" %in% class(override.ci.low) &&
               length(override.ci.low) == 1 &&
               override.ci.low == 0) {
      ci.low <- models[[i]]@ci.low
    } else if (is.numeric(override.ci.low) &&
               length(models) == 1 &&
               length(override.ci.low) == length(models[[i]]@coef)) {
      ci.low <- override.ci.low
    } else if (!"list" %in% class(override.ci.low)) {
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
    } else if (!is.numeric(override.ci.low[[i]])) {
      warning("Lower CIs provided for model", i, "are not numeric. Using default lower CIs.")
      ci.low <- models[[i]]@ci.low
    } else {
      ci.low <- override.ci.low[[i]]
    }
    models[[i]]@ci.low <- ci.low

    # upper bound of confidence intervals
    if (is.null(override.ci.up)) {
      # do nothing
    } else if (!"list" %in% class(override.ci.up) &&
               length(override.ci.up) == 1 &&
               override.ci.up == 0) {
      ci.up <- models[[i]]@ci.up
    } else if (is.numeric(override.ci.up) &&
               length(models) == 1 &&
               length(override.ci.up) == length(models[[i]]@coef)) {
      ci.up <- override.ci.up
    } else if (!"list" %in% class(override.ci.up)) {
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
    } else if (!is.numeric(override.ci.up[[i]])) {
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

#' Reorder a matrix vertically according to a vector of new positions
#'
#' Reorder a matrix vertically according to a vector of new positions.
#'
#' This function takes a matrix and reorders its rows based on a vector of new
#' positions.
#'
#' @param mat Input matrix.
#' @param new.order Vector of integer numbers with the new order of rows. The
#'   new order must contain as many elements as the matrix has rows and it must
#'   not contain NA values, duplicate entries, or gaps.
#' @return Reordered matrix.
#'
#' @keywords internal
#' @author Philip Leifeld
#' @seealso \code{\link{matrixreg}}
reorder <- function(mat, new.order) {
  if (is.null(new.order)) {
    return(mat)
  } else if (nrow(mat) != length(new.order)) {
    stop(paste("Error when reordering matrix: there are", nrow(mat),
               "rows, but you provided", length(new.order), "numbers."))
  } else if ("list" %in% class(new.order)) {
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

#' Constructor for \linkS4class{texreg} objects
#'
#' Constructor for \linkS4class{texreg} objects.
#'
#' This function creates a \linkS4class{texreg} object. A \linkS4class{texreg}
#' object contains information about coefficients, standard errors, p-values
#' (optional), and about goodness-of-fit statistics. Instead of standard
#' errors and p-values, a \linkS4class{texreg} object may also contain upper and
#' lower bounds of a confidence interval. \linkS4class{texreg} objects are used
#' by the \code{\link{texreg}} function to create LaTeX tables and other
#' representations of the model results.
#'
#' @param coef.names The names for the covariates in a model as a
#'   \code{character} vector (= row names).
#' @param coef The coefficients as a \code{numeric} vector. Can have length
#'   zero.
#' @param se The standard errors as a \code{numeric} vector. Can have length
#'   zero.
#' @param pvalues The p-values as a \code{numeric} vector. Can have length zero.
#' @param ci.low The lower bounds of the confidence intervals as a
#'   \code{numeric} vector. Can have length zero.
#' @param ci.up The upper bounds of the confidence intervals as a
#'   \code{numeric} vector. Can have length zero.
#' @param gof.names Names of the goodness-of-fit statistics as a
#'   \code{character} vector. Can have length zero.
#' @param gof Goodness-of-fit statistics as a \code{numeric} vector. Can have
#'   length zero.
#' @param gof.decimal A \code{logical} vector with as many elements as the
#'   \code{gof} argument, indicating whether the respective GOF statistic is a
#'   double (\code{TRUE}) or integer (\code{FALSE}) number.
#' @param model.name A name for the statistical model. Can be a \code{character}
#'   vector of length zero if there is no model name.
#' @return A \linkS4class{texreg} object representing the statistical model.
#'
#' @author Philip Leifeld
#' @seealso \code{\link{extract}}
#'
#' @references Leifeld, Philip (2013). texreg: Conversion of Statistical Model
#'   Output in R to LaTeX and HTML Tables. Journal of Statistical Software
#'   55(8): 1-24. \doi{10.18637/jss.v055.i08}.
#'
#' @examples
#' library("nlme")  # load library for fitting linear mixed effects models
#' model <- lme(distance ~ age, data = Orthodont, random = ~ 1)  # estimate
#' coefficient.names <- rownames(summary(model)$tTable)  # extract coef names
#' coefficients <- summary(model)$tTable[, 1]  # extract coefficient values
#' standard.errors <- summary(model)$tTable[, 2]  # extract standard errors
#' significance <- summary(model)$tTable[, 5]  #extract p-values
#'
#' lik <- summary(model)$logLik  # extract log likelihood
#' aic <- summary(model)$AIC  # extract AIC
#' bic <- summary(model)$BIC  # extract BIC
#' n <- nobs(model)  # extract number of observations
#' gof <- c(aic, bic, lik, n)  # create a vector of GOF statistics
#' gof.names <- c("AIC", "BIC", "Log Likelihood", "Num. obs.")  # names of GOFs
#' decimal.places <- c(TRUE, TRUE, TRUE, FALSE)  # last one is a count variable
#'
#' # create the texreg object
#' tr <- createTexreg(coef.names = coefficient.names,
#'                    coef = coefficients,
#'                    se = standard.errors,
#'                    pvalues = significance,
#'                    gof.names = gof.names,
#'                    gof = gof,
#'                    gof.decimal = decimal.places)
#'
#' @importFrom methods new
#' @export
createTexreg <- function(coef.names,
                         coef,
                         se = numeric(0),
                         pvalues = numeric(0),
                         ci.low = numeric(0),
                         ci.up = numeric(0),
                         gof.names = character(0),
                         gof = numeric(0),
                         gof.decimal = logical(0),
                         model.name = character(0)) {
  methods::new("texreg",
      coef.names = coef.names,
      coef = coef,
      se = se,
      pvalues = pvalues,
      ci.low = ci.low,
      ci.up = ci.up,
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal,
      model.name = model.name)
}

#' An S4 class to represent a statistical model as a texreg object
#'
#' An S4 class to represent a statistical model as a texreg object.
#'
#' A \linkS4class{texreg} object stores details about a statistical model. It
#' can be used for creating regression tables using \code{\link{screenreg}},
#' \code{\link{texreg}}, and similar functions.
#'
#' @slot coef.names The covariate names.
#' @slot coef The coefficients.
#' @slot se The standard errors.
#' @slot pvalues The p-values.
#' @slot ci.low The lower bounds of the confidence intervals.
#' @slot ci.up The upper bounds of the confidence intervals.
#' @slot gof.names The names of the goodness-of-fit statistics.
#' @slot gof The goodness-of-fit statistics.
#' @slot gof.decimal A vector describing for each GOF statistic whether it is a
#'   decimal value (\code{TRUE}) or an integer value (\code{FALSE}).
#' @slot model.name An optional model name. Can be of length zero.
#'
#' @author Philip Leifeld
#' @seealso \code{\link{extract}} \code{\link{createTexreg}}
#'
#' @references Leifeld, Philip (2013). texreg: Conversion of Statistical Model
#'   Output in R to LaTeX and HTML Tables. Journal of Statistical Software
#'   55(8): 1-24. \doi{10.18637/jss.v055.i08}.
#'
#' @export
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

#' Show method for pretty output of \linkS4class{texreg} objects
#'
#' Show method for pretty output of \linkS4class{texreg} objects.
#'
#' Print the different slots of \linkS4class{texreg} objects to the screen.
#'
#' @param object The \linkS4class{texreg} object to display.
#'
#' @author Philip Leifeld
#' @seealso \code{\link{extract}}, \code{\link{createTexreg}},
#'   \code{\link{screenreg}}
#'
#' @references Leifeld, Philip (2013). texreg: Conversion of Statistical Model
#'   Output in R to LaTeX and HTML Tables. Journal of Statistical Software
#'   55(8): 1-24. \doi{10.18637/jss.v055.i08}.
#'
#' @importFrom methods show
#' @export
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
