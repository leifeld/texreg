\name{texreg}
\alias{texreg}
\alias{htmlreg}
\alias{screenreg}
\alias{write.texreg}
\title{Convert regression output to LaTeX or HTML tables}
\description{Conversion of R regression output to LaTeX or HTML tables.}
\usage{
texreg(l, single.row=FALSE, no.margin=TRUE, leading.zero=TRUE, 
    table=TRUE, sideways=FALSE, float.pos="", stars=TRUE, 
    strong.signif=FALSE, symbol="\\\\cdot", use.packages=TRUE, 
    caption="Statistical models", label="table:coefficients", 
    dcolumn=TRUE, booktabs=TRUE, scriptsize=FALSE, 
    custom.names=NA, model.names=NA, digits=2, override.coef=0, 
    override.se=0, override.pval=0, omit.coef=NA, file=NA, 
    return.string=FALSE, ...)
htmlreg(l, single.row=FALSE, leading.zero=TRUE, stars=TRUE, 
    strong.signif=FALSE, symbol="&middot;", 
    caption="Statistical models", custom.names=NA, 
    model.names=NA, digits=2, override.coef=0, override.se=0, 
    override.pval=0, omit.coef=NA, file=NA, 
    return.string=FALSE, ...)
screenreg(l, single.row=FALSE, leading.zero=TRUE, stars=TRUE, 
    strong.signif=FALSE, custom.names=NA, model.names=NA, 
    digits=2, outer.rule="=", inner.rule="-", column.spacing=2, 
    override.coef=0, override.se=0, override.pval=0, 
    omit.coef=NA, file=NA, return.string=FALSE, ...)
}
\details{
texreg converts coefficients, standard errors, significance stars, 
and goodness-of-fit statistics of statistical models into LaTeX 
tables or HTML tables or into nicely formatted screen output for 
the R console. A list of several models can be combined in a 
single table. The output is customizable. New model types can be 
easily implemented.
}
\arguments{
\item{l}{ A statistical model or a list of statistical models. Lists of models can be specified as \code{l = list(model.1, model.2, ...)}. Different object types can also be mixed. }
\item{single.row}{ By default, a model parameter takes up two lines of the table: the standard error is listed in parentheses under the coefficient. This saves a lot of horizontal space on the page and is the default table format in most academic journals. If \code{single.row=TRUE} is activated, however, both coefficient and standard error are placed in a single table cell in the same line. }
\item{no.margin}{ In order to save space, inner margins of tables are switched off by default. To reactivate the default table spacing, set \code{no.margin=FALSE}. }
\item{leading.zero}{ Most journals require leading zeros of coefficients and standard errors (for example, \code{0.35}). This is also the default texreg behavior. Some journals, however, require omission of leading zeros (for example, \code{.35}). This can be achieved by setting \code{leading.zero=FALSE}. }
\item{table}{ By default, texreg puts the actual \code{tabular} object in a \code{table} floating environment. To get only the \code{tabular} object without the whole table header, set \code{table=FALSE}. }
\item{sideways}{ If \code{sideways=TRUE} is set, the \code{table} floating environment is replaced by a \code{sidewaystable} float, and the \code{rotating} package is loaded in the preamble. The argument only has an effect if \code{table=TRUE} is also set. }
\item{float.pos}{ This argument specifies where the table should be located on the page or in the document. By default, no floating position is specified, and LaTeX takes care of the position automatically. Possible values include \code{h} (here), \code{p} (page), \code{t} (top), \code{b} (bottom), any combination thereof, e.g. \code{tb}, or any of these values followed by an exclamation mark, e.g. \code{t!}, in order to enforce this position. The square brackets do not have to be specified. }
\item{stars}{ Should significance stars be printed after coefficients/standard errors, and should a note about significance levels be printed below the table? }
\item{strong.signif}{ Conventional regression tables report 99 percent significance as three stars, 95 percent as two stars, and 90 percent as one star. For ERGMs, however, the de facto standard is to report 99.9 percent as three stars, 99 percent as two stars, 95 percent as one star, and 90 percent as a single dot. By default, texreg follows conventional significance levels. To use the stronger ERGM standard, set \code{strong.signif=TRUE}. }
\item{symbol}{ By default, p values smaller than 0.1 and larger than 0.05 are reported as centered dots (\code{symbol="\\\\cdot"} for the LaTeX dot or \code{symbol="&middot;"} for the HTML dot) when \code{strong.signif=TRUE} is set. Any other mathematical LaTeX symbol or plain text symbol can be used, for example \code{symbol="\\\\circ"} for a small circle in LaTeX. }
\item{use.packages}{ If this argument is set to \code{TRUE} (= the default behavior), the required LaTeX packages are loaded in the beginning. If set to \code{FALSE}, the use package statements are omitted from the output. }
\item{caption}{ Set the caption of the table. }
\item{label}{ Set the label of the \code{table} environment. }
\item{dcolumn}{ Use the \code{dcolumn} LaTeX package to get a nice alignment of the coefficients. }
\item{booktabs}{ Use the \code{booktabs} LaTeX package to get thick horizontal rules in the output table. }
\item{scriptsize}{ To save horizontal space on the page, the table can be set in script size instead of normal text size by setting \code{scriptsize=TRUE}. }
\item{custom.names}{ By default, texreg uses the coefficient names which are stored in the models. The \code{custom.names} argument can be used to replace them by other character strings in the order of appearance. For example, if a table shows a total of three different independent variables (including the intercept), the argument \code{custom names=c("Intercept", "variable 1", "variable 2")} will replace their names in this order.

Sometimes it happens that the same variable has a different name in different models. In this case, the user can use this function to assign identical names. If possible, the rows will then be merged into a single row unless both rows contain values in the same column. }
\item{model.names}{ A string vector of labels for the models. By default, the models are named Model 1, Model 2, etc. Specifying \code{model.names=c("My name 1", "My name 2")} etc. overrides the default behavior. }
\item{digits}{ Set the number of decimal places for coefficients, standard errors and goodness-of-fit statistics. Do not use negative values! }
\item{override.coef}{ Set custom values for the coefficients. New coefficients are provided as a list of numeric vectors. The list contains vectors of coefficients for each model. There must be as many vectors of coefficients as there are models. For example, if there are two models with three model terms each, the argument could be specified as \code{override.coef=list(c(0.1, 0.2, 0.3), c(0.05, 0.06, 0.07))}. If there is only one model, custom values can be provided as a plain vector (not embedded in a list). For example: \code{override.coef=c(0.05, 0.06, 0.07)}. }
\item{override.se}{ Set custom values for the standard errors. New standard errors are provided as a list of numeric vectors. The list contains vectors of standard errors for each model. There must be as many vectors of standard errors as there are models. For example, if there are two models with three coefficients each, the argument could be specified as \code{override.se=list(c(0.1, 0.2, 0.3), c(0.05, 0.06, 0.07))}. If there is only one model, custom values can be provided as a plain vector (not embedded in a list). For example: \code{override.se=c(0.05, 0.06, 0.07)}. }
\item{override.pval}{ Set custom values for the p values. New p values are provided as a list of numeric vectors. The list contains vectors of p values for each model. There must be as many vectors of p values as there are models. For example, if there are two models with three coefficients each, the argument could be specified as \code{override.pval=list(c(0.1, 0.2, 0.3), c(0.05, 0.06, 0.07))}. If there is only one model, custom values can be provided as a plain vector (not embedded in a list). For example: \code{override.pval=c(0.05, 0.06, 0.07)}. }
\item{omit.coef}{ A character string which is used as a regular expression to remove coefficient rows from the table. For example, \code{omit.coef="group"} deletes all coefficient rows from the table where the name of the coefficient contains the character sequence "group". }
\item{file}{ Using this argument, the resulting table is written to a file rather than to the R prompt. The file name can be specified as a character string. Writing a table to a file can be useful for working with MS Office or LibreOffice. For example, an HTML table can be written to a file with the extension \code{.doc} and opened with MS Word. The table can then be simply copied into any Word document, retaining the formatting of the table. Note that LibreOffice can import only plain HTML; CSS decorations are not supported; the resulting tables do not retain the full formatting in LibreOffice. }
\item{return.string}{ Should the table be returned as a character object, which can later be printed again using the \code{cat()} function? }
\item{outer.rule}{ The character which is used to draw the outer horizontal line above and below a table (only for the \code{screenreg()} function). If an empty character object is provided (i.e., \code{outer.rule=""}), there will be no outer horizontal lines. Recommended values are \code{""}, \code{"="}, \code{"-"}, \code{"_"}, or \code{"#"}. }
\item{inner.rule}{ The character which is used to draw the inner horizontal line above and below a table (only for the \code{screenreg()} function). If an empty character object is provided (i.e., \code{outer.rule=""}), there will be no inner horizontal lines. Recommended values are \code{""}, \code{"-"}, or \code{"_"}. }
\item{column.spacing}{ The amount of space between any two columns of a table (only for the \code{screenreg()} function). By default, two spaces are used. If the tables do not fit on a single page horizontally, the value can be set to \code{1} or \code{0}. }
\item{...}{ Custom options to be passed on to the extract function. }
}
\seealso{
\code{\link{texreg-package} \link{extract} \link{extract-methods}}
}
\author{
Philip Leifeld (\url{http://www.philipleifeld.de})
}
\examples{
#Linear mixed-effects models
library(nlme)
model.1 <- lme(distance ~ age, data = Orthodont, random = ~ 1)
model.2 <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1)
table.string <- texreg(list(model.1, model.2))
cat(table.string)

#Ordinary least squares model (example from the 'lm' help file)
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2,10,20, labels=c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
table.string <- texreg(lm.D9)
cat(table.string)

#Create a 'fake' Office document containing a regression table
table.string <- htmlreg(list(model.1, model.2), file="texreg.doc")
unlink("texreg.doc")
}
\keyword{Programming|programming|interface}
\keyword{Programming|IO}
