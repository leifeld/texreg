# The texreg package was written by Philip Leifeld.
# Please use the issue tracker at http://github.com/leifeld/texreg
# for bug reports, help, or feature requests.

# matrixreg function

#' Convert regression output to LaTeX or HTML tables
#' 
#' Conversion of R regression output to LaTeX or HTML tables.
#'
#' @param l A statistical model or a list of statistical models. Lists of
#'  models can be specified as `l = list(model.1, model.2, ...)`.
#'  Different object types can also be mixed.
#' @param single.row By default, a model parameter takes up two lines of the
#'   table: the standard error is listed in parentheses under the coefficient.
#'   This saves a lot of horizontal space on the page and is the default table
#'   format in most academic journals. If `single.row = TRUE` is activated,
#'   however, both coefficient and standard error are placed in a single table
#'   cell in the same line.
#' @param stars The significance levels to be used to draw stars. Between 0 and
#'   4 threshold values can be provided as a numeric vector. For example,
#'   `stars = numeric(0)` will not print any stars and will not print any
#'   note about significance levels below the table. `stars = 0.05` will
#'   attach one single star to all coefficients where the p value is below 0.05.
#'   `stars = c(0.001, 0.01, 0.05, 0.1)` will print one, two, or three
#'   stars, or a symbol as specified by the `symbol` argument depending on
#'   the p values.
#' @param custom.model.names A character vector of labels for the models. By
#'   default, the models are named Model 1, Model 2, etc. Specifying
#'   `model.names = c("My name 1", "My name 2")` etc. overrides the default
#'   behavior.
#' @param custom.coef.names By default, \pkg{texreg} uses the coefficient names
#'   which are stored in the models. The `custom.coef.names` argument can
#'   be used to replace them by other character strings in the order of
#'   appearance. For example, if a table shows a total of three different
#'   coefficients (including the intercept), the argument
#'   `custom.coef.names = c("Intercept", "variable 1", "variable 2")` will
#'   replace their names in this order. 
#'   
#'   Sometimes it happens that the same
#'   variable has a different name in different models. In this case, the user
#'   can use this function to assign identical names. If possible, the rows will
#'   then be merged into a single row unless both rows contain values in the
#'   same column. 
#'   
#'   Where the argument contains an `NA` value, the original
#'   name of the coefficient is kept. For example, 
#'   `custom.coef.names = c(NA, "age", NA)` will only replace the second coef 
#'   name and leave the first and third name as they are in the original model.
#' @param custom.coef.map The `custom.coef.map` argument can be used to
#'   select, omit, rename, and reorder coefficients. Users must supply a named
#'   list of this form: 
#'   `list('x' = 'First variable', 'y' = NA, 'z' = Third variable')`. 
#'   With that particular example of `custom.coef.map`:
#'   
#'   1. coefficients will presented in order: x, y, z. 
#'   2. variable x will appear as "First variable", variable y will appear as 
#'     "y",and variable "z" will appear as "Third variable".
#'   3. all variables not named "x", "y", or "z" will be omitted from the table.
#' @param custom.gof.names A character vector which is used to replace the names
#'   of the goodness-of-fit statistics at the bottom of the table. The vector
#'   must have the same length as the number of GOF statistics in the final
#'   table. The argument works like the `custom.coef.names` argument, but
#'   for the GOF values. `NA` values can be included where the original GOF
#'   name should be kept.
#' @param custom.gof.rows A named list of vectors for new lines at the beginning
#'   of the GOF block of the table. For example, 
#'   `list("Random effects" = c("YES", "YES", "NO"), Observations = c(25, 25, 26))` 
#'   would insert two new
#'   rows into the table, at the beginning of the GOF block (i.e., after the
#'   coefficients). The rows can contain integer, numeric, or character objects.
#'   Note that this argument is processed after the `custom.gof.names`
#'   argument (meaning `custom.gof.names` should not include any of the new
#'   GOF rows) and before the `reorder.gof` argument (meaning that the new
#'   GOF order specified there should contain values for the new custom GOF
#'   rows). Arguments for custom columns are not affected because they only
#'   insert columns into the coefficient block.
#' @param digits Set the number of decimal places for coefficients, standard
#'   errors and goodness-of-fit statistics. Do not use negative values! The
#'   argument works like the `digits` argument in the `round` function
#'   of the \pkg{base} package.
#' @param leading.zero Most journals require leading zeros of coefficients and
#'   standard errors (for example, `0.35`). This is also the default texreg
#'   behavior. Some journals, however, require omission of leading zeros (for
#'   example, `.35`). This can be achieved by setting `leading.zero =
#'   FALSE`.
#' @param star.symbol 
#' @param symbol If four threshold values are handed over to the `stars`
#'   argument, p values smaller than the largest threshold value but larger than
#'   the second-largest threshold value are denoted by this symbol. The default
#'   symbol is `"\\\\cdot"` for the LaTeX dot, `"&middot;"` for the
#'   HTML dot, or simply `"."` for the ASCII dot. If the `texreg`
#'   function is used, any other mathematical LaTeX symbol or plain text symbol
#'   can be used, for example `symbol = "\\\\circ"` for a small circle
#'   (note that backslashes must be escaped). If the `htmlreg` function is
#'   used, any other HTML character or symbol can be used. For the
#'   `screenreg` function, only plain text characters can be used.
#' @param override.coef Set custom values for the coefficients. New coefficients
#'   are provided as a list of numeric vectors. The list contains vectors of
#'   coefficients for each model. There must be as many vectors of coefficients
#'   as there are models. For example, if there are two models with three model
#'   terms each, the argument could be specified as 
#'   `override.coef = list(c(0.1, 0.2, 0.3), c(0.05, 0.06, 0.07))`. If there is 
#'   only one model,
#'   custom values can be provided as a plain vector (not embedded in a list).
#'   For example: `override.coef = c(0.05, 0.06, 0.07)`.
#' @param override.se Set custom values for the standard errors. New standard
#'   errors are provided as a list of numeric vectors. The list contains vectors
#'   of standard errors for each model. There must be as many vectors of
#'   standard errors as there are models. For example, if there are two models
#'   with three coefficients each, the argument could be specified as
#'   `override.se = list(c(0.1, 0.2, 0.3), c(0.05, 0.06, 0.07))`. If there
#'   is only one model, custom values can be provided as a plain vector (not
#'   embedded in a list). For example: `override.se = c(0.05, 0.06, 0.07)`.
#'   Overriding standard errors can be useful for the implementation of robust
#'   SEs, for example.
#' @param override.pvalues Set custom values for the p values. New p values are
#'   provided as a list of numeric vectors. The list contains vectors of p
#'   values for each model. There must be as many vectors of p values as there
#'   are models. For example, if there are two models with three coefficients
#'   each, the argument could be specified as `override.pvalues =
#'   list(c(0.1, 0.2, 0.3), c(0.05, 0.06, 0.07))`. If there is only one model,
#'   custom values can be provided as a plain vector (not embedded in a list).
#'   For example: `override.pvalues = c(0.05, 0.06, 0.07)`. Overriding p
#'   values can be useful for the implementation of robust SEs and p values, for
#'   example.
#' @param override.ci.low Set custom lower confidence interval bounds. This
#'   works like the other override arguments, with one exception: if confidence
#'   intervals are provided here and in the `override.ci.up` argument, the
#'   standard errors and p values as well as the `ci.force` argument are
#'   ignored.
#' @param override.ci.up Set custom upper confidence interval bounds. This works
#'   like the other override arguments, with one exception: if confidence
#'   intervals are provided here and in the `override.ci.low` argument, the
#'   standard errors and p values as well as the `ci.force` argument are
#'   ignored.
#' @param omit.coef A character string which is used as a regular expression to
#'   remove coefficient rows from the table. For example, `omit.coef =
#'   "group"` deletes all coefficient rows from the table where the name of the
#'   coefficient contains the character sequence "group". More complex regular
#'   expressions can be used to filter out several kinds of model terms, for
#'   example `omit.coef = "(thresh)|(ranef)"` to remove all model terms
#'   matching either "thresh" or "ranef".The `omit.coef` argument is
#'   processed after the `custom.coef.names` argument, so the regular
#'   expression should refer to the custom coefficient names. To omit GOF
#'   entries instead of coefficient entries, use the custom arguments of the
#'   extract functions instead (see the help entry of the \link{extract}
#'   function or \link{extract-methods}.
#' @param reorder.coef Reorder the rows of the coefficient block of the
#'   resulting table in a custom way. The argument takes a vector of the same
#'   length as the number of coefficients. For example, if there are three
#'   coefficients, `reorder.coef = c(3, 2, 1)` will put the third
#'   coefficient in the first row and the first coefficient in the third row.
#'   Reordering can be sensible because interaction effects are often added to
#'   the end of the model output although they were specified earlier in the
#'   model formula. Note: Reordering takes place after processing custom
#'   coefficient names and after omitting coefficients, so the
#'   `custom.coef.names` and `omit.coef` arguments should follow the
#'   original order.
#' @param reorder.gof Reorder the rows of the goodness-of-fit block of the
#'   resulting table in a custom way. The argument takes a vector of the same
#'   length as the number of GOF statistics. For example, if there are three
#'   goodness-of-fit rows, `reorder.gof = c(3, 2, 1)` will exchange the
#'   first and the third row. Note: Reordering takes place after processing
#'   custom GOF names and after adding new custom GOF rows, so the
#'   `custom.gof.names` and `custom.gof.rows` arguments should follow
#'   the original order, and the `reorder.gof` argument should contain
#'   values for any rows that are added through the `custom.gof.rows`
#'   argument.
#' @param ci.force Should confidence intervals be used instead of the default
#'   standard errors and p values? Most models implemented in the \pkg{texreg}
#'   package report standard errors and p values by default while few models
#'   report confidence intervals. However, the functions in the \pkg{texreg}
#'   package can convert standard errors and into confidence intervals if
#'   desired. To enforce confidence intervals instead of standard errors, the
#'   `ci.force` argument accepts either a logical value indicating whether
#'   all models or none of the models should be forced to report confidence
#'   intervals (`ci.force = TRUE` for all and `ci.force = FALSE` for
#'   none) or a vector of logical values indicating for each model separately
#'   whether the model should be forced to report confidence intervals (e.g.,
#'   `ci.force = c(FALSE, TRUE, FALSE)`). Confidence intervals are computed
#'   using the standard normal distribution (z values based on the `qnorm`
#'   function).
#' @param ci.force.level If the `ci.force` argument is used to convert
#'   standard errors to confidence intervals, what confidence level should be
#'   used? By default, `0.95` is used (i.e., an alpha value of 0.05).
#' @param ci.test If confidence intervals are reported, the `ci.test`
#'   argument specifies the reference value to establish whether a
#'   coefficient/CI is significant. The default value `ci.test = 0`, for
#'   example, will attach a significance star to coefficients if the confidence
#'   interval does not contain `0`. If no star should be printed at all,
#'   `ci.test = NULL` can be used. The `ci.test` argument works both
#'   for models with native support for confidence intervals and in cases where
#'   the `ci.force` argument is used.
#' @param bold [only in the `texreg` and `htmlreg` functions] The p
#'   value threshold below which the coefficient shall be formatted in a bold
#'   font. For example, `bold = 0.05` will cause all coefficients which are
#'   significant at the 95\% level to be formatted in bold. Note that this is
#'   not compatible with the dcolumn argument in the `texreg` function. If
#'   both are `TRUE`, dcolumn is switched off and a warning message
#'   appears. Note also that it is advisable to use `stars = FALSE`
#'   together with the `bold` argument because having both bolded
#'   coefficients and significance stars usually does not make any sense.
#' @param groups This argument can be used to group the rows of the table into
#'   blocks. For example, there could be one block for hypotheses and another
#'   block for control variables. Each group has a heading, and the row labels
#'   within a group are indented. The partitions must be handed over as a list
#'   of named numeric vectors, where each number is a row index and each name is
#'   the heading of the group. Example: `groups = list("first group" = 1:4,
#'   "second group" = 7:8)`.
#' @param custom.columns An optional list of additional text columns to be
#'   inserted into the table, for example coefficient types. The list should
#'   contain one or more character vectors with as many character or numeric
#'   elements as there are rows. If the vectors in the list are named, the names
#'   are used as labels in the table header. For example, `custom.columns =
#'   list(type = c("a", "b", "c"), 1:3)` will add two columns; the first one is
#'   labeled while the second one is not. Note that the numeric elements of the
#'   second column will be converted to character objects in this example. The
#'   consequence is that decimal alignment with the \pkg{dcolumn} package is
#'   switched off in these columns. Note that this argument is processed after
#'   any arguments that affect the number of rows.
#' @param custom.col.pos An optional integer vector of positions for the columns
#'   given in the `custom columns` argument. For example, if there are
#'   three custom columns, `custom.col.pos = c(1, 3, 3)` will insert the
#'   first custom column before the first column of the original table and the
#'   remaining two custom columns after the second column of the original table.
#'   By default, all custom columns are placed after the first column, which
#'   usually contains the coefficient names.
#' @param dcolumn only in the `texreg` function] Use the `dcolumn` 
#'   LaTeX package to get a nice alignment of the coefficients (recommended).
#' @param ... Custom options to be passed on to the extract function. For
#'   example, most extract methods provide custom options for the inclusion or
#'   exclusion of specific goodness-of-fit statistics. See the help entries of
#'   [extract()] and [extract-methods] for more information.
#'  
#' @details
#' texreg converts coefficients, standard errors, significance stars, 
#' and goodness-of-fit statistics of statistical models into LaTeX 
#' tables or HTML tables or into nicely formatted screen output for 
#' the R console. A list of several models can be combined in a 
#' single table. The output is customizable. New model types can be 
#' easily implemented. Confidence intervals can be used instead of 
#' standard errors and p-values.
#' 
#' The `texreg()` function creates LaTeX code for inclusion 
#' in a LaTeX document or for usage with \pkg{Sweave} or \pkg{knitr}.
#' 
#' The `htmlreg()` function creates HTML code. Tables in HTML 
#' The `htmlreg()` function creates HTML code. Tables in HTML 
#' format can be saved with a ".html" extension and displayed in 
#' a web browser. Alternatively, they can be saved with a ".doc" 
#' extension and opened in MS Word for inclusion in office 
#' documents. `htmlreg()`` also works with \pkg{knitr} and HTML 
#' or Markdown. Note that the `inline.css`, `doctype`, 
#' `html.tag`, `head.tag`, and `body.tag` arguments 
#' must be adjusted for the different purposes (see the description 
#'   of the arguments).
#' 
#' The `screenreg()` function creates text representations of 
#' tables and prints them to the R console. This is an alternative 
#' to the `summary` method and serves easy model comparison. 
#' Moreover, once a table has been prepared in the R console, it 
#' can be later exported to LaTeX or HTML with little extra effort 
#' because the majority of arguments of the three functions is 
#' identical.
#' 
#' The `wordreg()` function creates a Microsoft Word document with the
#' requested table. 
#' 
#' The `matrixreg()` function creates text representations of 
#' tables and prints them to the R console. 
#' 
#' The `csvreg()` function creates a spreadsheet in CSV containing the
#' information in the tables. 
#' 
#' The `huxtablereg()` function creates a 
#' [huxtable::huxtable()] object using the 'huxtable'
#' package. This allows output to HTML, LaTeX, Word, Excel, Powerpoint, and 
#' RTF. The object can be formatted using huxtable package functions. See also
#' [huxtable::huxreg()].
#' 
#' @export
#'
#' @examples
#' # Linear mixed-effects models
#' library(nlme)
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
#' # Create a 'fake' Office document containing a regression table
#' htmlreg(list(model.1, model.2), file = "texreg.doc", 
#'   inline.css = FALSE, doctype = TRUE, html.tag = TRUE, 
#'   head.tag = TRUE, body.tag = TRUE)
#' unlink("texreg.doc")
#' 
#' if (requireNamespace('huxtable')) {
#'   hr <- huxtablereg(list(model.1, model.2))
#'   hr <- huxtable::set_bottom_border(hr, 1, -1, 0.4)
#'   hr <- huxtable::set_bold(hr, 1:nrow(hr), 1, TRUE)
#'   hr <- huxtable::set_bold(hr, 1, -1, TRUE)
#'   hr <- huxtable::set_all_borders(hr, 4, 2, 0.4)
#'   hr <- huxtable::set_all_border_colors(hr, 4, 2, "red")
#'   hr
#'   \dontrun{
#'     huxtable::quick_pdf(hr)
#'     huxtable::quick_docx(hr)
#'     # or use in a knitr document
#'   }
#' }
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
                      ...) {
  
  # unnamed arguments to environment
  dots <- list(...)
  
  # argument for internal use
  if (!"output.type" %in% names(dots)) {
    dots[["output.type"]] <- "ascii"
  }
  
  # extract
  models <- get.data(l, ...)  # extract relevant coefficients, SEs, GOFs etc.
  models <- override(models, override.coef, override.se, override.pvalues, 
                     override.ci.low, override.ci.up)
  if (dots$output.type != "latex") {
    models <- tex.replace(models, type = "screen")  # convert TeX code to text
  }
  models <- ciforce(models, ci.force = ci.force, ci.level = ci.force.level)
  gof.names <- get.gof(models)  # extract names of GOFs (before adding any GOF rows)
  models <- correctDuplicateCoefNames(models)
  
  # arrange coefficients and GOFs nicely in a matrix
  if (dots$output.type == "latex") {
    gof.matrix <- aggregate.matrix(models = models,
                                   gof.names = gof.names,
                                   custom.gof.names = custom.gof.names,
                                   custom.gof.rows = custom.gof.rows,
                                   reorder.gof = reorder.gof,
                                   digits = digits,
                                   leading.zero = leading.zero,
                                   latex = TRUE,
                                   dcolumn = dcolumn,
                                   returnobject = "gof.matrix")
    m <- aggregate.matrix(models = models,
                          gof.names = gof.names,
                          custom.gof.names = custom.gof.names,
                          custom.gof.rows = custom.gof.rows,
                          reorder.gof = reorder.gof,
                          digits = digits,
                          leading.zero = leading.zero,
                          latex = TRUE,
                          dcolumn = dcolumn,
                          returnobject = "m")
  } else {
    gof.matrix <- aggregate.matrix(models = models,
                                   gof.names = gof.names,
                                   custom.gof.names = custom.gof.names,
                                   custom.gof.rows = custom.gof.rows,
                                   reorder.gof = reorder.gof,
                                   digits = digits,
                                   leading.zero = leading.zero,
                                   latex = FALSE,
                                   dcolumn = dcolumn,
                                   returnobject = "gof.matrix")
    m <- aggregate.matrix(models = models,
                          gof.names = gof.names,
                          custom.gof.names = custom.gof.names,
                          custom.gof.rows = custom.gof.rows,
                          reorder.gof = reorder.gof,
                          digits = digits,
                          leading.zero = leading.zero,
                          latex = FALSE,
                          dcolumn = dcolumn,
                          returnobject = "m")
  }
  
  if (!is.null(custom.coef.map)) {
    m <- custommap(m, custom.coef.map)
  } else {
    m <- omit_rename(m,
                     omit.coef = omit.coef,
                     custom.coef.names = custom.coef.names)
  }
  m <- rearrangeMatrix(m)  # resort matrix and conflate duplicate entries
  m <- as.data.frame(m)
  
  mod.names <- modelnames(l, models, custom.model.names)  # model names
  
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
  
  # output matrix
  if (dots$output.type == "ascii") {
    output.matrix <- outputmatrix(m, single.row, neginfstring = "-Inf", 
                                  posinfstring = "Inf", leading.zero, digits, 
                                  se.prefix = " (", se.suffix = ")", star.prefix = " ", star.suffix = "", 
                                  stars, dcolumn = TRUE, star.symbol = star.symbol, symbol = symbol, 
                                  bold = bold, bold.prefix = "", bold.suffix = "", ci = ci, 
                                  ci.test = ci.test)
  } else if (dots$output.type == "latex") {
    output.matrix <- outputmatrix(m, single.row, 
                                  neginfstring = "\\multicolumn{1}{c}{$-\\infty$}", 
                                  posinfstring = "\\multicolumn{1}{c}{$\\infty$}", leading.zero, digits, 
                                  se.prefix = " \\; (", se.suffix = ")", star.prefix = "^{", 
                                  star.suffix = "}", stars, dcolumn = dcolumn, star.symbol = star.symbol,
                                  symbol = symbol, bold = bold, bold.prefix = "\\mathbf{", 
                                  bold.suffix = "}", ci = ci, semicolon = ";\\ ", ci.test = ci.test, rowLabelType = 'latex')
  } else if (dots$output.type == "html") {
    output.matrix <- outputmatrix(m, single.row, neginfstring = "-Inf", 
                                  posinfstring = "Inf", leading.zero, digits, 
                                  se.prefix = " (", se.suffix = ")", star.symbol = star.symbol, 
                                  star.prefix = paste0("<sup", dots$css.sup, ">"), 
                                  star.suffix = "</sup>", stars, dcolumn = TRUE, symbol = symbol, 
                                  bold = bold, bold.prefix = "<b>", bold.suffix = "</b>", ci = ci, 
                                  ci.test = ci.test)
  }
  
  # grouping
  if (dots$output.type == "latex") {
    output.matrix <- grouping(output.matrix, groups, indentation = "    ", 
                              single.row = single.row, prefix = "", suffix = "", rowLabelType = "latex")
  } else {
    output.matrix <- grouping(output.matrix, groups, indentation = "    ", 
                              single.row = single.row, prefix = "", suffix = "", rowLabelType = "text")
  }
  
  # combine the coefficient and gof matrices vertically
  coef.names <- output.matrix[output.matrix[, 1] != "", 1]
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
  if (dots$output.type == "ascii") {
    output.matrix <- rbind(c("", mod.names), output.matrix)
  }
  
  # add custom columns
  output.matrix <- customcolumns(output.matrix, custom.columns, custom.col.pos, 
                                 single.row = single.row, numcoef = nrow(m), groups = groups, 
                                 modelnames = TRUE)
  
  # attributes required for printing functions
  if ("include.attributes" %in% names(dots)) {
    if (dots$include.attributes) {
      attr(output.matrix, "ci") <- ci
      attr(output.matrix, "ci.test") <- ci.test
      attr(output.matrix, "gof.names") <- gof.matrix[, 1]
      attr(output.matrix, "coef.names") <- coef.names
      attr(output.matrix, "mod.names") <- mod.names
    }
  }
  
  return(output.matrix)
} 

# screenreg function
#' @rdname matrixreg
#' @param file Using this argument, the resulting table is written to a file rather than to
#' the R prompt. The file name can be specified as a character string. Writing a table to a
#' file can be useful for working with MS Office or LibreOffice. For example, using the
#' `htmlreg()` function, an HTML table can be written to a file with the extension 
#' `.doc` and opened with MS Word. The table can then be simply copied into any Word
#' document, retaining the formatting of the table. Note that LibreOffice can import only
#' plain HTML; CSS decorations are not supported; the resulting tables do not retain the
#' full formatting in LibreOffice.
#' @export
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
                             ...)
  
  gof.names <- attr(output.matrix, 'gof.names')
  coef.names <- attr(output.matrix, 'coef.names')
  mod.names <- attr(output.matrix, 'mod.names')
  ci <- attr(output.matrix, 'ci')
  ci.test <- attr(output.matrix, 'ci.test')
  
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
#' @rdname matrixreg
#' @export
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
                             include.attributes = TRUE,
                             output.type = "latex",
                             ...)
  
  gof.names <- attr(output.matrix, "gof.names")
  coef.names <- attr(output.matrix, "coef.names")
  mod.names <- attr(output.matrix, "mod.names")
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
#' @rdname matrixreg
#' @export
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
#' @rdname matrixreg
#' @export
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
                   caption = 'Statistical Models',
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
#' @rdname matrixreg
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

# huxtablereg function
#' @rdname matrixreg
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
  mx <- eval(mr.call)
  
  mx <- trimws(mx)
  
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
