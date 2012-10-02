\name{texreg}
\alias{texreg}
\title{Convert regression output to LaTeX tables}
\description{Conversion of R regression output to LaTeX tables.}
\usage{
texreg(l, single.row=FALSE, no.margin=TRUE, leading.zero=TRUE, 
    table=TRUE, sideways=FALSE, float.pos="", 
    strong.signif=FALSE, symbol="\\\\cdot", use.packages=TRUE, 
    caption="Statistical models", label="table:coefficients", 
    dcolumn=TRUE, booktabs=TRUE, scriptsize=FALSE, 
    custom.names=NA, model.names=NA, digits=2, ...)
}
\details{
texreg converts coefficients, standard errors, significance stars, 
and goodness-of-fit statistics of statistical models into LaTeX 
tables. Either a single model or a list of models can be combined 
in a single table. The output is customizable. Models can be handed 
over as lists of objects. New model types can be easily implemented.
}
\arguments{
\item{l}{ An lm, lme, lrm, gls, glm, ergm, plm, pmg, rq or clogit object. The method accepts a single model or a list of models. Lists of models can be specified as \code{l = list(model.1, model.2, ...)}. Different object types can also be mixed. }
\item{single.row}{ By default, a model parameter takes up two lines of the table: the standard error is listed in parentheses under the coefficient. This saves a lot of horizontal space on the page and is the default table format in most academic journals. If \code{single.row=TRUE} is activated, however, both coefficient and standard error are placed in a single table cell in the same line. }
\item{no.margin}{ In order to save space, inner margins of tables are switched off by default. To reactivate the default table spacing, set \code{no.margin=FALSE}. }
\item{leading.zero}{ Most journals require leading zeros of coefficients and standard errors (for example, \code{0.35}). This is also the default texreg behavior. Some journals, however, require omission of leading zeros (for example, \code{.35}). This can be achieved by setting \code{leading.zero=FALSE}. }
\item{table}{ By default, texreg puts the actual \code{tabular} object in a \code{table} floating environment. To get only the \code{tabular} object without the whole table header, set \code{table=FALSE}. }
\item{sideways}{ If \code{sideways=TRUE} is set, the \code{table} floating environment is replaced by a \code{sidewaystable} float, and the \code{rotating} package is loaded in the preamble. The argument only has an effect if \code{table=TRUE} is also set. }
\item{float.pos}{ This argument specifies where the table should be located on the page or in the document. By default, no floating position is specified, and LaTeX takes care of the position automatically. Possible values include \code{h} (here), \code{p} (page), \code{t} (top), \code{b} (bottom), any combination thereof, e.g. \code{tb}, or any of these values followed by an exclamation mark, e.g. \code{t!}, in order to enforce this position. The square brackets do not have to be specified. }
\item{strong.signif}{ Conventional regression tables report 99 percent significance as three stars, 95 percent as two stars, and 90 percent as one star. For ERGMs, however, the de facto standard is to report 99.9 percent as three stars, 99 percent as two stars, 95 percent as one star, and 90 percent as a single dot. By default, texreg follows conventional significance levels. To use the stronger ERGM standard, set \code{strong.signif=TRUE}. }
\item{symbol}{ By default, p values smaller than 0.1 and larger than 0.05 are reported as dots (\code{symbol="\\\\cdot"}) when \code{strong.signif=TRUE} is set. Any other mathematical LaTeX symbol can be used, for example \code{symbol="\\\\circ"}. }
\item{use.packages}{ If this argument is set to \code{TRUE} (= the default behavior), the required LaTeX packages are loaded in the beginning. If set to \code{FALSE}, the use package statements are omitted from the output. }
\item{caption}{ Set the caption of the \code{table} environment. }
\item{label}{ Set the label of the \code{table} environment. }
\item{dcolumn}{ Use the \code{dcolumn} LaTeX package to get a nice alignment of the coefficients. }
\item{booktabs}{ Use the \code{booktabs} LaTeX package to get thick horizontal rules in the output table. }
\item{scriptsize}{ To save horizontal space on the page, the table can be set in script size instead of normal text size by setting \code{scriptsize=TRUE}. }
\item{custom.names}{ By default, texreg uses the coefficient names which are stored in the models. The \code{custom.names} argument can be used to replace them by other character strings in the order of appearance. For example, if a table shows a total of three different independent variables (including the intercept), the argument \code{custom names=c("Intercept", "variable 1", "variable 2")} will replace their names in this order.

Sometimes it happens that the same variable has a different name in different models. In this case, the user can use this function to assign identical names. If possible, the rows will then be merged into a single row unless both rows contain values in the same column. }
\item{model.names}{ A string vector of labels for the models. By default, the models are named Model 1, Model 2, etc. Specifying \code{model.names=c("My name 1", "My name 2")} etc. overrides the default behavior. }
\item{digits}{ Set the number of decimal places for coefficients, standard errors and goodness-of-fit statistics. Do not use negative values! }
\item{...}{ Custom options to be passed on to the extract function. }
}
\seealso{
\code{\link{texreg-package}}
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

#Exponential random graph models
library(network)
library(ergm)
mat <- rbinom(400,1,0.16)  #create a matrix
mat <- matrix(mat, nrow=20)
nw <- network(mat)
model.1 <- ergm(nw ~ edges)
model.2 <- ergm(nw ~ edges + mutual)
model.3 <- ergm(nw ~ edges + mutual + twopath)
table.string <- texreg(list(model.1, model.2, model.3))
cat(table.string)

#Ordinary least squares model (example from the 'lm' help file)
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2,10,20, labels=c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
table.string <- texreg(lm.D9)
cat(table.string)
}
\keyword{Programming|programming|interface}
\keyword{Programming|IO}
