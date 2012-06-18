\name{texreg}
\alias{texreg}
\title{Get network data}
\description{Conversion of R regression output to LaTeX tables.}
\usage{
texreg(l, single.row=FALSE, no.margin=TRUE, leading.zero=TRUE, 
    table=TRUE, strong.signif=TRUE, symbol="\\\\cdot", 
    use.packages=TRUE, caption="Statistical models", 
    label="table:coefficients", dcolumn=TRUE, booktabs=TRUE, 
    scriptsize=TRUE)
}
\details{
texreg converts coefficients, standard errors, significance stars, 
and goodness-of-fit statistics of lme objects (that is, linear mixed 
effects models) and ergm objects (that is, exponential random graph models) 
into LaTeX tables. Either a single model or a list of models can be 
combined in a single table. The output is customizable. Models can be handed 
over as lists of objects.
}
\arguments{
\item{l}{ An ergm or lme object, which was fitted using the ergm() command in the ergm package or using the lme() command in the nlme package. The method accepts a single model or a list of models. Lists of models can be specified as \code{l = list(model.1, model.2, ...)}. ergm and lme objects can also be mixed. }
\item{single.row}{ By default, a model parameter takes up two lines of the table: the standard error is listed in parentheses under the coefficient. This saves a lot of horizontal space on the page and is the default table format in most academic journals. If \code{single.row=TRUE} is activated, however, both coefficient and standard error are placed in a single table cell in the same line. }
\item{no.margin}{ In order to save space, inner margins of tables are switched off by default. To reactivate the default table spacing, set \code{no.margin=FALSE}. }
\item{leading.zero}{ Most journals require leading zeros of coefficients and standard errors (for example, \code{0.35}). This is also the default texreg behavior. Some journals, however, require omission of leading zeros (for example, \code{.35}). This can be achieved by setting \code{leading.zero=FALSE}. }
\item{table}{ By default, texreg puts the actual \code{tabular} object in a \code{table} floating environment. To get only the \code{tabular} object without the whole table header, set \code{table=FALSE}. }
\item{strong.signif}{ Conventional regression tables report 99 percent significance as three stars, 95 percent as two stars, and 90 percent as one star. For ERGMs, however, the de facto standard is to report 99.9 percent as three stars, 99 percent as two stars, 95 percent as one star, and 90 percent as a single dot. By default, texreg follows the stronger ERGM standard. To obtain conventional significance levels, set \code{strong.signif=FALSE}. }
\item{symbol}{ By default, p values smaller than 0.1 and larger than 0.05 are reported as dots (\code{symbol="\\\\cdot"}). Any other mathematical LaTeX symbol can be used, for example \code{symbol="\\\\circ"}. }
\item{use.packages}{ If this argument is set to \code{TRUE} (= the default behavior), the required LaTeX packages are loaded in the beginning. If set to \code{FALSE}, the use package statements are omitted from the output. }
\item{caption}{ Set the caption of the \code{table} environment. }
\item{label}{ Set the label of the \code{table} environment. }
\item{dcolumn}{ Use the \code{dcolumn} LaTeX package to get a nice alignment of the coefficients. }
\item{booktabs}{ Use the \code{booktabs} LaTeX package to get thick horizontal rules in the output table. }
\item{scriptsize}{ To save horizontal space on the page, the table is set in script size instead of normal text size by default. If \code{scriptsize=FALSE}, the default font size is used. }
}
\seealso{
\code{\link{texreg-package}}
}
\author{
Philip Leifeld (\url{http://www.philipleifeld.de})
}
\examples{
library(nlme)
model.1 <- lme(distance ~ age, data = Orthodont, random = ~ 1)
model.2 <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1)
table.string <- texreg(list(model.1, model.2))
cat(table.string)
}
\keyword{Programming|programming|interface}
\keyword{Programming|IO}
