\name{plotreg}
\alias{plotreg}
\alias{coefplot}
\title{Create coefficient plots from statistical model output}
\description{Create coefficient plots from statistical model output.}
\usage{

plotreg(l, file = NULL, custom.model.names = NULL, 
    custom.coef.names = NULL, custom.note = NULL, 
    override.coef = 0, override.se = 0, override.pval = 0, 
    override.ci.low = 0, override.ci.up = 0, 
    omit.coef = NULL, reorder.coef = NULL, ci.level = 0.95, 
    use.se = FALSE, mfrow = TRUE, xlim = NULL, cex = 2.5, 
    lwd.zerobar = 4, lwd.vbars = 1, lwd.inner = 7, 
    lwd.outer = 5, ylab.cex = 1.0, signif.light = "#fbc9b9", 
    signif.medium = "#f7523a", signif.dark = "#bd0017", 
    insignif.light = "#c5dbe9", insignif.medium = "#5a9ecc", 
    insignif.dark = "#1c5ba6", ...)

coefplot(labels, estimates, lower.inner = NULL, 
    upper.inner = NULL, lower.outer = NULL, 
    upper.outer = NULL, signif.outer = TRUE, 
    xlab = "Coefficients and confidence intervals", 
    main = "Coefficient plot", xlim = NULL, 
    cex = 2.5, lwd.zerobar = 4, lwd.vbars = 1, 
    lwd.inner = 7, lwd.outer = 5, ylab.cex = 1.0, 
    signif.light = "#fbc9b9", signif.medium = "#f7523a", 
    signif.dark = "#bd0017", insignif.light = "#c5dbe9", 
    insignif.medium = "#5a9ecc", insignif.dark = "#1c5ba6", 
    ...)
}
\details{
The \code{coefplot} function produces coefficient plots (i.e., forest plots applied to point estimates and confidence intervals). It accepts raw data (the lower and upper bounds of inner and outer confidence intervals as well as the point estimates and their names) as input data. Significant coefficients and intervals can be plotted in a different color.

The \code{plotreg} function is a wrapper for the \code{coefplot} function and works much like the 
\code{\link{screenreg}}, \code{\link{texreg}}, and \code{\link{htmlreg}} functions. It accepts a single or multiple statistical models as input and internally extracts the relevant data from the models. If confidence intervals are not defined in the extract method of a statistical model (see \link{extract} and \link{extract-methods}), the default standard errors are converted to confidence intervals. Most of the arguments work either like in the \code{\link{screenreg}}, \code{\link{texreg}}, and \code{\link{htmlreg}} functions, or they work like in the \code{coefplot} function.
}
\arguments{
\item{l}{A statistical model or a list of statistical models. Lists of models can be specified as \code{l = list(model.1, model.2, ...)}. Different object types can also be mixed.}
\item{file}{Using this argument, the resulting table is written to a file rather than to the R prompt. The file name can be specified as a character string. The file extension is automatically recognized. \code{pdf}, \code{ps}, \code{png}, \code{bmp}, \code{jpg}, and \code{tiff} are supported.}
\item{custom.model.names}{A character vector of labels for the models. By default, the models are named Model 1, Model 2, etc. Specifying \code{model.names = c("My name 1", "My name 2")} etc. overrides the default behavior.}
\item{custom.coef.names}{By default, \code{plotreg} uses the coefficient names which are stored in the models. The \code{custom.coef.names} argument can be used to replace them by other character strings in the order of appearance. For example, if a model shows a total of three coefficients (including the intercept), the argument \code{custom.coef.names = c("Intercept", "variable 1", "variable 2")} will replace their names in this order.}
\item{custom.note}{With this argument, a replacement text for the \code{xlab} note below the diagram can be provided. If an empty character object is provided (\code{custom.note = ""}), the note will be omitted completely.}
\item{override.coef}{Set custom values for the coefficients. New coefficients are provided as a list of numeric vectors. The list contains vectors of coefficients for each model. There must be as many vectors of coefficients as there are models. For example, if there are two models with three model terms each, the argument could be specified as \code{override.coef = list(c(0.1, 0.2, 0.3), c(0.05, 0.06, 0.07))}. If there is only one model, custom values can be provided as a plain vector (not embedded in a list). For example: \code{override.coef = c(0.05, 0.06, 0.07)}.}
\item{override.se}{Set custom values for the standard errors. This only has an effect where standard errors are converted into confidence intervals because no other CIs are present. New standard errors are provided as a list of numeric vectors. The list contains vectors of standard errors for each model. There must be as many vectors of standard errors as there are models. For example, if there are two models with three coefficients each, the argument could be specified as \code{override.se = list(c(0.1, 0.2, 0.3), c(0.05, 0.06, 0.07))}. If there is only one model, custom values can be provided as a plain vector (not embedded in a list). For example: \code{override.se = c(0.05, 0.06, 0.07)}. Overriding standard errors can be useful for the implementation of robust SEs, for example.}
\item{override.pval}{Set custom values for the p values. This only has an effect where standard errors are converted into confidence intervals because no other CIs are present. In this case, significance is derived from the p values rather than the confidence intervals. New p values are provided as a list of numeric vectors. The list contains vectors of p values for each model. There must be as many vectors of p values as there are models. For example, if there are two models with three coefficients each, the argument could be specified as \code{override.pval = list(c(0.1, 0.2, 0.3), c(0.05, 0.06, 0.07))}. If there is only one model, custom values can be provided as a plain vector (not embedded in a list). For example: \code{override.pval = c(0.05, 0.06, 0.07)}. Overriding p values can be useful for the implementation of robust SEs and p values, for example.}
\item{override.ci.low}{Set custom lower confidence interval bounds. This works like the other override arguments, with one exception: if confidence intervals are provided here and in the \code{override.ci.up} argument, the standard errors and p values as well as the \code{ci.force} argument are ignored.}
\item{override.ci.up}{Set custom upper confidence interval bounds. This works like the other override arguments, with one exception: if confidence intervals are provided here and in the \code{override.ci.low} argument, the standard errors and p values as well as the \code{ci.force} argument are ignored.}
\item{omit.coef}{A character string which is used as a regular expression to remove coefficient rows from the table. For example, \code{omit.coef = "group"} deletes all coefficient rows from the diagram where the name of the coefficient contains the character sequence "group". More complex regular expressions can be used to filter out several kinds of model terms, for example \code{omit.coef = "(thresh)|(ranef)"} to remove all model terms matching either "thresh" or "ranef".The \code{omit.coef} argument is processed after the \code{custom.coef.names} argument, so the regular expression should refer to the custom coefficient names.}
\item{reorder.coef}{Reorder the rows of the coefficient block of the resulting table in a custom way. The argument takes a vector of the same length as the number of coefficients. For example, if there are three coefficients, \code{reorder.coef = c(3, 2, 1)} will put the third coefficient in the first row and the first coefficient in the third row. Reordering can be sensible because interaction effects are often added to the end of the model output although they were specified earlier in the model formula. Note: Reordering takes place after processing custom coefficient names and after omitting coefficients, so the \code{custom.coef.names} and \code{omit.coef} arguments should follow the original order.}
\item{ci.level}{If standard errors are converted to confidence intervals (because a model does not natively support CIs), what confidence level should be used for the outer confidence interval? By default, \code{0.95} is used (i.e., an alpha value of 0.05).}
\item{use.se}{Use one standard error for the inner horizontal bar and two standard errors from the estimate for the outer horizontal bar (instead of confidence intervals). Only available if standard errors can be extracted from the model using the respective \code{\link{extract}} function.}
\item{mfrow}{If multiple models are handed over as the \code{l} argument, several plots are produced. If \code{mfrow = TRUE} is set, multiple diagrams are aligned on the same page. If \code{mfrow = FALSE} is set, each diagram per model comes out as a separate plot.}
\item{xlim}{Horizontal limits. In the \code{coefplot} function, they must be provided as a vector with two numeric, e.g., \code{xlim = c(-5, 5)} for displaying a range from \code{-5} to \code{+5}. In the \code{plotreg} function, they can be provided either as such a vector with two values or as a list of vectors (with each entry corresponding to a model in \code{l}). }
\item{lwd.zerobar}{ Line width of the vertical gray bar at the \code{x} value of \code{0}. To remove the line, set \code{lwd.zerobar = 0}. }
\item{lwd.vbars}{ Line width of the thin vertical gray bars. To remove them completely, set \code{lwd.vbars = 0}. }

\item{labels}{The names of the model terms. They are arranged on the left axis.}
\item{estimates}{The coefficients (point estimates) of the model terms. They are depicted as bold dots in the coefficient plot.}
\item{lower.inner}{The lower bounds of the inner confidence intervals, provided as a vector. Inner CI means more relaxed (lower confidence level, higher alpha) because fewer observations have to fall into the CI, therefore the CI gets smaller.}
\item{upper.inner}{The upper bounds of the inner confidence intervals, provided as a vector. Inner CI means more relaxed (lower confidence level, higher alpha) because fewer observations have to fall into the CI, therefore the CI gets smaller.}
\item{lower.outer}{The lower bounds of the outer confidence intervals, provided as a vector. Outer CI means stricter or narrower (higher confidence level, lower alpha) because more observations have to fall into the CI, therefore the CI gets larger.}
\item{upper.outer}{The upper bounds of the outer confidence intervals, provided as a vector. Outer CI means stricter or narrower (higher confidence level, lower alpha) because more observations have to fall into the CI, therefore the CI gets larger.}

\item{signif.outer}{Different colors are used for significant estimates and confidence intervals. If \code{signif.outer = TRUE}, the outer CIs are used to evaluate significance, otherwise the inner CIs are used.}
\item{xlab}{The label of the x axis.}
\item{main}{The main title or heading of the plot.}
\item{cex}{Size of the point representing the estimate.}
\item{lwd.inner}{Line width of the inner confidence interval.}
\item{lwd.outer}{Line width of the outer confidence interval.}
\item{ylab.cex}{Size of the coefficient labels. The size of the x axis labels can be adjusted by using argument \code{cex.axis = 1.0}.}
\item{signif.light}{Color of outer confidence intervals for significant model terms.}
\item{signif.medium}{Color of inner confidence intervals for significant model terms.}
\item{signif.dark}{Color of point estimates and labels for significant model terms.}
\item{insignif.light}{Color of outer confidence intervals for insignificant model terms.}
\item{insignif.medium}{Color of inner confidence intervals for insignificant model terms.}
\item{insignif.dark}{Color of point estimates and labels for insignificant model terms.}

\item{...}{Custom options to be passed on to the extract function or the graphics device. See the help entries of \link{extract} and \link{extract-methods} for more information.}
}
\seealso{
\code{\link{texreg-package} \link{extract} \link{extract-methods} \link{texreg}}
}
\author{
Philip Leifeld (\url{http://www.philipleifeld.com})
}
\examples{
#example from the 'lm' help file:
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2,10,20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
screenreg(lm.D9)  # print model output to the R console
plotreg(lm.D9)    # plot model output as a diagram
}
\keyword{print}
\keyword{misc}
\keyword{utilities}
\keyword{IO}
\keyword{programming|interface}
