\name{texreg}
\alias{texreg}
\alias{htmlreg}
\alias{screenreg}
\title{Convert regression output to LaTeX or HTML tables}
\description{Conversion of R regression output to LaTeX or HTML tables.}
\usage{
texreg(l, file = NULL, single.row = FALSE, stars = c(0.001, 
    0.01, 0.05), custom.model.names = NULL, 
    custom.coef.names = NULL, custom.gof.names = NULL, 
    custom.note = NULL, digits = 2, leading.zero = TRUE, 
    symbol = "\\\\cdot", override.coef = 0, override.se = 0, 
    override.pval = 0, override.ci.low = 0, override.ci.up = 0, 
    omit.coef = NULL, reorder.coef = NULL, reorder.gof = NULL, 
    ci.force = FALSE, ci.force.level = 0.95, ci.test = 0, 
    groups = NULL, bold = 0.00, center = TRUE, 
    caption = "Statistical models", caption.above = FALSE, 
    label = "table:coefficients", booktabs = FALSE, 
    dcolumn = FALSE, sideways = FALSE, longtable = FALSE, 
    use.packages = TRUE, table = TRUE, no.margin = TRUE, 
    fontsize = NULL, scalebox = NULL, float.pos = "", ...)

htmlreg(l, file = NULL, single.row = FALSE, stars = c(0.001, 
    0.01, 0.05), custom.model.names = NULL, 
    custom.coef.names = NULL, custom.gof.names = NULL, 
    custom.note = NULL, digits = 2, leading.zero = TRUE, 
    symbol = "&middot;", override.coef = 0, override.se = 0, 
    override.pval = 0, override.ci.low = 0, override.ci.up = 0, 
    omit.coef = NULL, reorder.coef = NULL, reorder.gof = NULL, 
    ci.force = FALSE, ci.force.level = 0.95, ci.test = 0, 
    groups = NULL, bold = 0.00, center = TRUE, 
    caption = "Statistical models", caption.above = FALSE, 
    star.symbol = "*", inline.css = TRUE, doctype = TRUE, 
    html.tag = FALSE, head.tag = FALSE, body.tag = FALSE, 
    indentation = "", vertical.align.px = 0, ...)

screenreg(l, file = NULL, single.row = FALSE, stars = c(0.001, 
    0.01, 0.05), custom.model.names = NULL, 
    custom.coef.names = NULL, custom.gof.names = NULL, 
    custom.note = NULL, digits = 2, leading.zero = TRUE, 
    symbol = ".", override.coef = 0, override.se = 0, 
    override.pval = 0, override.ci.low = 0, override.ci.up = 0, 
    omit.coef = NULL, reorder.coef = NULL, reorder.gof = NULL, 
    ci.force = FALSE, ci.force.level = 0.95, ci.test = 0, 
    groups = NULL, column.spacing = 2, outer.rule = "=", 
    inner.rule = "-", ...)
}
\details{
texreg converts coefficients, standard errors, significance stars, 
and goodness-of-fit statistics of statistical models into LaTeX 
tables or HTML tables or into nicely formatted screen output for 
the R console. A list of several models can be combined in a 
single table. The output is customizable. New model types can be 
easily implemented. Confidence intervals can be used instead of 
standard errors and p values.

The \code{texreg()} function creates LaTeX code for inclusion 
in a LaTeX document or for usage with \pkg{Sweave} or \pkg{knitr}.

The \code{htmlreg()} function creates HTML code. Tables in HTML 
format can be saved with a ".html" extension and displayed in 
a web browser. Alternatively, they can be saved with a ".doc" 
extension and opened in MS Word for inclusion in office 
documents. \code{htmlreg()} also works with \pkg{knitr} and HTML 
or Markdown. Note that the \code{inline.css}, \code{doctype}, 
\code{html.tag}, \code{head.tag}, and \code{body.tag} arguments 
must be adjusted for the different purposes (see the description 
of the arguments).

The \code{screenreg()} function creates text representations of 
tables and prints them to the R console. This is an alternative 
to the \code{summary} method and serves easy model comparison. 
Moreover, once a table has been prepared in the R console, it 
can be later exported to LaTeX or HTML with little extra effort 
because the majority of arguments of the three functions is 
identical.
}
\arguments{
\item{l}{A statistical model or a list of statistical models. Lists of models can be specified as \code{l = list(model.1, model.2, ...)}. Different object types can also be mixed.}
\item{file}{Using this argument, the resulting table is written to a file rather than to the R prompt. The file name can be specified as a character string. Writing a table to a file can be useful for working with MS Office or LibreOffice. For example, using the \code{htmlreg} function, an HTML table can be written to a file with the extension \code{.doc} and opened with MS Word. The table can then be simply copied into any Word document, retaining the formatting of the table. Note that LibreOffice can import only plain HTML; CSS decorations are not supported; the resulting tables do not retain the full formatting in LibreOffice.}
\item{single.row}{By default, a model parameter takes up two lines of the table: the standard error is listed in parentheses under the coefficient. This saves a lot of horizontal space on the page and is the default table format in most academic journals. If \code{single.row = TRUE} is activated, however, both coefficient and standard error are placed in a single table cell in the same line.}
\item{stars}{The significance levels to be used to draw stars. Between 0 and 4 threshold values can be provided as a numeric vector. For example, \code{stars = numeric(0)} will not print any stars and will not print any note about significance levels below the table. \code{stars = 0.05} will attach one single star to all coefficients where the p value is below 0.05. \code{stars = c(0.001, 0.01, 0.05, 0.1)} will print one, two, or three stars, or a symbol as specified by the \code{symbol} argument depending on the p values.}
\item{custom.model.names}{A character vector of labels for the models. By default, the models are named Model 1, Model 2, etc. Specifying \code{model.names = c("My name 1", "My name 2")} etc. overrides the default behavior.}
\item{custom.coef.names}{By default, \pkg{texreg} uses the coefficient names which are stored in the models. The \code{custom.coef.names} argument can be used to replace them by other character strings in the order of appearance. For example, if a table shows a total of three different coefficients (including the intercept), the argument \code{custom.coef.names = c("Intercept", "variable 1", "variable 2")} will replace their names in this order.

Sometimes it happens that the same variable has a different name in different models. In this case, the user can use this function to assign identical names. If possible, the rows will then be merged into a single row unless both rows contain values in the same column.

Where the argument contains an \code{NA} value, the original name of the coefficient is kept. For example, \code{custom.coef.names = c(NA, "age", NA)} will only replace the second coef name and leave the first and third name as they are in the original model.}
\item{custom.gof.names}{A character vector which is used to replace the names of the goodness-of-fit statistics at the bottom of the table. The vector must have the same length as the number of GOF statistics in the final table. The argument works like the \code{custom.coef.names} argument, but for the GOF values. \code{NA} values can be included where the original GOF name should be kept.}
\item{custom.note}{With this argument, a replacement text for the significance note below the table can be provided. If an empty character object is provided (\code{custom.note = ""}), the note will be omitted completely. If some character string is provided (e.g., \code{custom.note = "My note"}), the significance legend is replaced by \code{My note}. The original significance legend can be included by inserting the \code{\%stars} wildcard. For example, a custom note can be added right after the significance legend by providing \code{custom.note = "\%stars. My note"}.}
\item{digits}{Set the number of decimal places for coefficients, standard errors and goodness-of-fit statistics. Do not use negative values! The argument works like the \code{digits} argument in the \code{round} function of the \pkg{base} package.}
\item{leading.zero}{Most journals require leading zeros of coefficients and standard errors (for example, \code{0.35}). This is also the default texreg behavior. Some journals, however, require omission of leading zeros (for example, \code{.35}). This can be achieved by setting \code{leading.zero = FALSE}.}
\item{symbol}{If four threshold values are handed over to the \code{stars} argument, p values smaller than the largest threshold value but larger than the second-largest threshold value are denoted by this symbol. The default symbol is \code{"\\\\cdot"} for the LaTeX dot, \code{"&middot;"} for the HTML dot, or simply \code{"."} for the ASCII dot. If the \code{texreg} function is used, any other mathematical LaTeX symbol or plain text symbol can be used, for example \code{symbol = "\\\\circ"} for a small circle (note that backslashes must be escaped). If the \code{htmlreg} function is used, any other HTML character or symbol can be used. For the \code{screenreg} function, only plain text characters can be used.}
\item{override.coef}{Set custom values for the coefficients. New coefficients are provided as a list of numeric vectors. The list contains vectors of coefficients for each model. There must be as many vectors of coefficients as there are models. For example, if there are two models with three model terms each, the argument could be specified as \code{override.coef = list(c(0.1, 0.2, 0.3), c(0.05, 0.06, 0.07))}. If there is only one model, custom values can be provided as a plain vector (not embedded in a list). For example: \code{override.coef = c(0.05, 0.06, 0.07)}.}
\item{override.se}{Set custom values for the standard errors. New standard errors are provided as a list of numeric vectors. The list contains vectors of standard errors for each model. There must be as many vectors of standard errors as there are models. For example, if there are two models with three coefficients each, the argument could be specified as \code{override.se = list(c(0.1, 0.2, 0.3), c(0.05, 0.06, 0.07))}. If there is only one model, custom values can be provided as a plain vector (not embedded in a list). For example: \code{override.se = c(0.05, 0.06, 0.07)}. Overriding standard errors can be useful for the implementation of robust SEs, for example.}
\item{override.pval}{Set custom values for the p values. New p values are provided as a list of numeric vectors. The list contains vectors of p values for each model. There must be as many vectors of p values as there are models. For example, if there are two models with three coefficients each, the argument could be specified as \code{override.pval = list(c(0.1, 0.2, 0.3), c(0.05, 0.06, 0.07))}. If there is only one model, custom values can be provided as a plain vector (not embedded in a list). For example: \code{override.pval = c(0.05, 0.06, 0.07)}. Overriding p values can be useful for the implementation of robust SEs and p values, for example.}
\item{override.ci.low}{Set custom lower confidence interval bounds. This works like the other override arguments, with one exception: if confidence intervals are provided here and in the \code{override.ci.up} argument, the standard errors and p values as well as the \code{ci.force} argument are ignored.}
\item{override.ci.up}{Set custom upper confidence interval bounds. This works like the other override arguments, with one exception: if confidence intervals are provided here and in the \code{override.ci.low} argument, the standard errors and p values as well as the \code{ci.force} argument are ignored.}
\item{omit.coef}{A character string which is used as a regular expression to remove coefficient rows from the table. For example, \code{omit.coef = "group"} deletes all coefficient rows from the table where the name of the coefficient contains the character sequence "group". More complex regular expressions can be used to filter out several kinds of model terms, for example \code{omit.coef = "(thresh)|(ranef)"} to remove all model terms matching either "thresh" or "ranef".The \code{omit.coef} argument is processed after the \code{custom.coef.names} argument, so the regular expression should refer to the custom coefficient names. To omit GOF entries instead of coefficient entries, use the custom arguments of the extract functions instead (see the help entry of the \link{extract} function or \link{extract-methods}.}
\item{reorder.coef}{Reorder the rows of the coefficient block of the resulting table in a custom way. The argument takes a vector of the same length as the number of coefficients. For example, if there are three coefficients, \code{reorder.coef = c(3, 2, 1)} will put the third coefficient in the first row and the first coefficient in the third row. Reordering can be sensible because interaction effects are often added to the end of the model output although they were specified earlier in the model formula. Note: Reordering takes place after processing custom coefficient names and after omitting coefficients, so the \code{custom.coef.names} and \code{omit.coef} arguments should follow the original order.}
\item{reorder.gof}{Reorder the rows of the goodness-of-fit block of the resulting table in a custom way. The argument takes a vector of the same length as the number of GOF statistics. For example, if there are three goodness-of-fit rows, \code{reorder.gof = c(3, 2, 1)} will exchange the first and the third row. Note: Reordering takes place after processing custom GOF names, so the \code{custom.gof.names} argument should follow the original order.}
\item{ci.force}{Should confidence intervals be used instead of the default standard errors and p values? Most models implemented in the \pkg{texreg} package report standard errors and p values by default while few models report confidence intervals. However, the functions in the \pkg{texreg} package can convert standard errors and into confidence intervals if desired. To enforce confidence intervals instead of standard errors, the \code{ci.force} argument accepts either a logical value indicating whether all models or none of the models should be forced to report confidence intervals (\code{ci.force = TRUE} for all and \code{ci.force = FALSE} for none) or a vector of logical values indicating for each model separately whether the model should be forced to report confidence intervals (e.g., \code{ci.force = c(FALSE, TRUE, FALSE)}).}
\item{ci.force.level}{If the \code{ci.force} argument is used to convert standard errors to confidence intervals, what confidence level should be used? By default, \code{0.95} is used (i.e., an alpha value of 0.05).}
\item{ci.test}{If confidence intervals are reported, the \code{ci.test} argument specifies the reference value to establish whether a coefficient/CI is significant. The default value \code{ci.test = 0}, for example, will attach a significance star to coefficients if the confidence interval does not contain \code{0}. If no star should be printed at all, \code{ci.test = NULL} can be used. The \code{ci.test} argument works both for models with native support for confidence intervals and in cases where the \code{ci.force} argument is used.}
\item{groups}{This argument can be used to group the rows of the table into blocks. For example, there could be one block for hypotheses and another block for control variables. Each group has a heading, and the row labels within a group are indented. The partitions must be handed over as a list of named numeric vectors, where each number is a row index and each name is the heading of the group. Example: \code{groups = list("first group" = 1:4, "second group" = 7:8)}.}

\item{bold}{[only in the \code{texreg} and \code{htmlreg} functions] The p value threshold below which the coefficient shall be formatted in a bold font. For example, \code{bold = 0.05} will cause all coefficients which are significant at the 95\% level to be formatted in bold. Note that this is not compatible with the dcolumn argument in the \code{texreg} function. If both are \code{TRUE}, dcolumn is switched off and a warning message appears. Note also that it is advisable to use \code{stars = FALSE} together with the \code{bold} argument because having both bolded coefficients and significance stars usually does not make any sense.}
\item{center}{[only in the \code{texreg} and \code{htmlreg} functions] Should the table be horizontally aligned at the center of the page?}
\item{caption}{[only in the \code{texreg} and \code{htmlreg} functions] Set the caption of the table. }
\item{caption.above}{[only in the \code{texreg} and \code{htmlreg} functions] Should the caption of the table be placed above the table? By default, it is placed below the table.}

\item{label}{[only in the \code{texreg} function] Set the label of the \code{table} environment.}
\item{booktabs}{[only in the \code{texreg} function] Use the \code{booktabs} LaTeX package to get thick horizontal rules in the output table (recommended).}
\item{dcolumn}{[only in the \code{texreg} function] Use the \code{dcolumn} LaTeX package to get a nice alignment of the coefficients (recommended).}
\item{sideways}{[only in the \code{texreg} function] If \code{sideways = TRUE} is set, the \code{table} floating environment is replaced by a \code{sidewaystable} float, and the \code{rotating} package is loaded in the preamble. The argument only has an effect if \code{table = TRUE} is also set.}
\item{longtable}{ [only in the \code{texreg} function] If \code{longtable = TRUE} is set, the \code{longtable} environment from the \code{longtable} LaTeX package is used to set tables across multiple pages. Note that this argument is not compatible with the \code{sideways} and \code{scalebox} arguments. These arguments will be automatically switched off when \code{longtable = TRUE} is set. }
\item{use.packages}{[only in the \code{texreg} function] If this argument is set to \code{TRUE} (= the default behavior), the required LaTeX packages are loaded in the beginning. If set to \code{FALSE}, the use package statements are omitted from the output.}
\item{table}{[only in the \code{texreg} function] By default, texreg puts the actual \code{tabular} object in a \code{table} floating environment. To get only the \code{tabular} object without the whole table header, set \code{table = FALSE}.}
\item{no.margin}{[only in the \code{texreg} function] In order to save space, inner margins of tables are switched off by default. To reactivate the default table spacing, set \code{no.margin = FALSE}. }
\item{fontsize}{[only in the \code{texreg} function] The \code{fontsize} argument serves to change the font size used in the table. Valid values are \code{"tiny"}, \code{"scriptsize"}, \code{"footnotesize"}, \code{"small"}, \code{"normalsize"}, \code{"large"}, \code{"Large"}, \code{"LARGE"}, \code{"huge"}, and \code{"Huge"}. Note that the \code{scalebox} argument often achieves better results when the goal is to change the size of the table. }
\item{scalebox}{[only in the \code{texreg} function] The \code{scalebox} argument serves to resize the table. For example, \code{scalebox = 1.0} is equivalent to the normal size, \code{scalebox = 0.5} decreases the size of the table by one half, and \code{scalebox = 2.0} doubles the space occupied by the table. Note that the \code{scalebox} argument does not work when the \code{longtable} argument is used. }
\item{float.pos}{[only in the \code{texreg} function] This argument specifies where the table should be located on the page or in the document. By default, no floating position is specified, and LaTeX takes care of the position automatically. Possible values include \code{h} (here), \code{p} (page), \code{t} (top), \code{b} (bottom), any combination thereof, e.g. \code{tb}, or any of these values followed by an exclamation mark, e.g. \code{t!}, in order to enforce this position. The square brackets do not have to be specified.}

\item{star.symbol}{[only in the \code{htmlreg} function] Alternative characters for the significance stars can be specified. This is useful if \pkg{knitr} and Markdown are used for HTML report generation. In Markdown, asterisks or stars are interpreted as special characters, so they have to be escaped. To make \code{htmlreg} compatible with Markdown, specify \code{star.symbol = "\\*"}. Note that some other modifications are recommended for usage with \pkg{knitr} in combination with Markdown or HTML (see the \code{inline.css}, \code{doctype}, \code{html.tag}, \code{head.tag}, and \code{body.tag} arguments).}
\item{inline.css}{[only in the \code{htmlreg} function] Should the CSS stylesheets be embedded directly in the code of the table (\code{inline.css = TRUE}), or should the CSS stylesheets be enclosed in the <head> tag, that is, separated from the table code (\code{inline.css = FALSE})? Having inline CSS code makes the code of the table more complex, but sometimes it may be helpful when only the table shall be printed, without the head of the HTML file (for example when the table is embedded in a \pkg{knitr} report). As a rule of thumb: use inline CSS if the table is not saved to a file.}
\item{doctype}{[only in the \code{htmlreg} function] Should the first line of the HTML code contain the DOCTYPE definition? If \code{TRUE}, the HTML 4 TRANSITIONAL version is used. If \code{FALSE}, no DOCTYPE will be included. Omitting the DOCTYPE can be helpful when the \pkg{knitr} package is used to generate HTML code because \pkg{knitr} requires only the plain table, not the whole HTML document including the document type declaration. Including the DOCTYPE can be helpful when the code is saved to a file, for example as an MS Word document.}
\item{html.tag}{[only in the \code{htmlreg} function] Should the table code (and possibly the <body> and <head> tags) be enclosed in an <html> tag? Suppressing this tag is recommended when \pkg{knitr} is used for dynamic HTML or Markdown report generation. Including this tag is recommended when the code is saved to a file, for example as an MS Word document.}
\item{head.tag}{[only in the \code{htmlreg} function] Should the <head> tag (including CSS definitions and title/caption) be included in the HTML code? Suppressing this tag is recommended when \pkg{knitr} is used for dynamic HTML or Markdown report generation. Including this tag is recommended when the code is saved to a file, for example as an MS Word document.}
\item{body.tag}{[only in the \code{htmlreg} function] Should the table code be enclosed in a <body> HTML tag? Suppressing this tag is recommended when \pkg{knitr} is used for dynamic HTML or Markdown report generation. Including this tag is recommended when the code is saved to a file, for example as an MS Word document.}
\item{indentation}{[only in the \code{htmlreg} function] Characters used for indentation of the HTML code. By default, \code{indentation = ""} uses no indentation. Any number of spaces or characters can be used instead. For example, \code{indentation = "  "} uses two spaces of (additional) indentation for each subelement.}
\item{vertical.align.px}{[only in the \code{htmlreg} function] Vertical alignment of significance stars. Browsers differ in their ways of displaying superscripted significance stars; in some browsers the stars are elevated by default, and in other browsers the stars are aligned vertically with the text, without any actual superscripting. This argument controls by how many additional pixels the stars are elevated. The default setting of 0 uses the defaults of the browser. In RStudio's internal browser, this looks OK, but in Firefox, this looks too low. A value of 4 looks OK in Firefox, for example, but is above the line in RStudio's internal browser. }

\item{column.spacing}{[only in the \code{screenreg} function] The amount of space between any two columns of a table. By default, two spaces are used. If the tables do not fit on a single page horizontally, the value can be set to \code{1} or \code{0}.}
\item{outer.rule}{[only in the \code{screenreg} function] The character which is used to draw the outer horizontal line above and below a table. If an empty character object is provided (i.e., \code{outer.rule = ""}), there will be no outer horizontal lines. Recommended values are \code{""}, \code{"="}, \code{"-"}, \code{"_"}, or \code{"#"}.}
\item{inner.rule}{[only in the \code{screenreg} function] The character which is used to draw the inner horizontal line above and below a table. If an empty character object is provided (i.e., \code{outer.rule = ""}), there will be no inner horizontal lines. Recommended values are \code{""}, \code{"-"}, or \code{"_"}.}

\item{...}{Custom options to be passed on to the extract function. For example, most extract methods provide custom options for the inclusion or exclusion of specific goodness-of-fit statistics. See the help entries of \link{extract} and \link{extract-methods} for more information.}
}
\references{
Leifeld, Philip (2013). texreg: Conversion of Statistical Model Output in R to 
LaTeX and HTML Tables. Journal of Statistical Software, 55(8), 1-24.
\url{http://www.jstatsoft.org/v55/i08/}.
}
\seealso{
\code{\link{texreg-package} \link{extract} \link{extract-methods} \link{plotreg}}
}
\author{
Philip Leifeld (\url{http://www.philipleifeld.de})
}
\examples{
#Linear mixed-effects models
library(nlme)
model.1 <- lme(distance ~ age, data = Orthodont, random = ~ 1)
model.2 <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1)
texreg(list(model.1, model.2), booktabs = TRUE, dcolumn = TRUE)

#Ordinary least squares model (example from the 'lm' help file)
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2,10,20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
table.string <- texreg(lm.D9, return.string = TRUE)
cat(table.string)

#Create a 'fake' Office document containing a regression table
htmlreg(list(model.1, model.2), file = "texreg.doc", 
    inline.css = FALSE, doctype = TRUE, html.tag = TRUE, 
    head.tag = TRUE, body.tag = TRUE)
unlink("texreg.doc")
}
\keyword{print}
\keyword{misc}
\keyword{utilities}
\keyword{IO}
\keyword{programming|interface}
