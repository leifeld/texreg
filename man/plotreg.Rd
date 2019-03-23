\name{plotreg}
\alias{plotreg}
\title{Create coefficient plots from statistical model output using ggplot2}
\description{Create coefficient plots from statistical model output using ggplot2.}
\usage{
plotreg(l, file = NULL, custom.model.names = NULL, custom.title = NULL,
        custom.coef.names = NULL, custom.note = NULL, override.coef = 0,
        override.se = 0, override.pval = 0, override.ci.low = 0, override.ci.up = 0,
        omit.coef = NULL, reorder.coef = NULL, ci.level = 0.95, use.se = FALSE,
        type = "facet", theme = theme_bw(), custom.coef.map = NULL, 
        signif.light = "#FBC9B9",signif.medium = "#F7523A", signif.dark = "#BD0017", 
        insignif.light = "#C5DBE9",insignif.medium = "#5A9ECC", 
        insignif.dark = "#1C5BA6", ggsave.output = NULL, ...)
}
\details{
The \code{plotreg} function produces coefficient plots (i.e., forest plots applied to point estimates and confidence intervals) and works much like the \code{\link{screenreg}}, \code{\link{texreg}}, \code{\link{htmlreg}}, \code{\link{matrixreg}}, \code{\link{csvreg}} and \code{\link{wordreg}} functions. It accepts a single or multiple statistical models as input and internally extracts the relevant data from the models. If confidence intervals are not defined in the extract method of a statistical model (see \link{extract} and \link{extract-methods}), the default standard errors are converted to confidence intervals. Most of the arguments work either like in the \code{\link{screenreg}}, \code{\link{texreg}}, and \code{\link{htmlreg}} \code{\link{matrixreg}}, \code{\link{csvreg}} and \code{\link{wordreg}} functions.It is possible to display the plots in two ways: using the \code{\link{facet}} option, one forest plot applied to point estimates and confidence intervals will be visualised in case there is only one model. If there is more than one model each one will be plotted next to the other; using the \code{\link{forest}} option, coefficients from one or more models will be grouped together and displayed as a forest plot. 
}
\arguments{
\item{l}{A statistical model or a list of statistical models. Lists of models can be specified as \code{l = list(model.1, model.2, ...)}. Different object types can also be mixed.}
\item{file}{Using this argument, the resulting table is written to a file rather than to the R prompt. The file name can be specified as a character string. The file extension is automatically recognized. \code{pdf}, \code{ps}, \code{png}, \code{bmp}, \code{jpg}, and \code{tiff} are supported.}
\item{custom.model.names}{A character vector of labels for the models. By default, the models are named Model 1, Model 2, etc. Specifying \code{model.names = c("My name 1", "My name 2")} etc. overrides the default behavior.}
\item{custom.title}{With this argument, a replacement text for the \code{ggtitle} that provides a title above the diagram can be provided. If an empty character object is provided (\code{custom.note = ""}), the note will be omitted completely.}
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
\item{type}{The default option is \code{\link{facet}}. If one model only is specified, it will print one forest plot applied to point estimates and confidence intervals. If more than one model is specified, it will print as many plots as the number of models in a column of plots. Alternatively, if \code{\link{forest}} is specified, coefficients from one or more models will be grouped together and displayed as a forest plot.}
\item{theme}{The \code{\link{theme}} argument can be used to customomise the non data component of your plot. The default \code{\link{theme}} is \code{theme_bw} and it can be substituted by any other theme compatible with \pkg{ggplot2}.}
\item{custom.coef.map}{The \code{custom.coef.map} argument can be used to select, omit, rename, and reorder coefficients. Users must supply a named list of this form: \code{list('x' = 'First variable', 'y' = NA, 'z' = 'Third variable')}. With that particular example of \code{custom.coef.map}, 
1. coefficients will presented in order: x, y, z.
2. variable x will appear as "First variable", variable y will appear as "y", and variable "z" will appear as "Third variable".
3. all variables not named "x", "y", or "z" will be omitted from the table.}
\item{signif.light}{Color of outer confidence intervals for significant model terms.}
\item{signif.medium}{Color of inner confidence intervals for significant model terms.}
\item{signif.dark}{Color of point estimates and labels for significant model terms.}
\item{insignif.light}{Color of outer confidence intervals for insignificant model terms.}
\item{insignif.medium}{Color of inner confidence intervals for insignificant model terms.}
\item{insignif.dark}{Color of point estimates and labels for insignificant model terms.}
\item{ggsave.output}{can be specified to directly save the plot in the desired format. E.g. \code{ggsave.output = "outputest.png"}}
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
plotreg(lm.D9, theme = theme_dark(),  ggtitle ="my title", ggsave.output = "myplot.pdf")
#plot model output as a diagram customising theme and title and automatically saving in pdf.
}
\keyword{print}
\keyword{misc}
\keyword{utilities}
\keyword{IO}
\keyword{programming|interface}
