#' Create coefficient plots from statistical model output using \code{\link{ggplot2}}.
#'
#' Create coefficient plots of \R regression output using \code{\link{ggplot2}}.
#'
#'  The \code{plotreg} function produces coefficient plots (i.e., forest plots 
#'  applied to point estimates and confidence intervals) and works much like 
#'  the \code{\link{screenreg}}, \code{\link{texreg}}, \code{\link{htmlreg}}, 
#'  \code{\link{matrixreg}} and \code{\link{wordreg}} functions. It accepts a 
#'  single or multiple statistical models as input and internally extracts the
#'  relevant data from the models. If confidence intervals are not defined 
#'  in the extract method of a statistical model (see \link{extract} and 
#'  \link{extract-methods}), the default standard errors are converted to 
#'  confidence intervals. Most of the arguments work either like in the
#'  \code{\link{screenreg}}, \code{\link{texreg}}, and \code{\link{htmlreg}} 
#'  \code{\link{matrixreg}} and \code{\link{wordreg}} functions.It is 
#'  possible to display the plots in two ways: using the \code{\link{facet}} 
#'  option, one forest plot applied to point estimates and confidence intervals 
#'  will be visualised in case there is only one model. If there is more than 
#'  one model each one will be plotted next to the other; using the 
#'  \code{\link{forest}} option, coefficients from one or more models will 
#'  be grouped together and displayed as a forest plot. 
#' 
#'
#' @param l A statistical model or a list of statistical models. Lists of
#'   models can be specified as \code{l = list(model.1, model.2, ...)}.
#'   Different object types can also be mixed.
#' @param custom.model.names A character vector of labels for the models. By
#'   default, the models are named "Model 1", "Model 2", etc. Specifying
#'   \code{model.names = c("My name 1", "My name 2")} etc. overrides the default
#'   behavior.
#' @param custom.title With this argument, a replacement text for the \code{ggtitle} 
#'   that provides a title above the diagram can be provided. If an empty character 
#'   object is provided (\code{custom.note = ""}), the note will be omitted completely.
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
#' @param custom.note With this argument, a replacement text for the \code{xlab} 
#'   note below the diagram can be provided. If an empty character object is provided
#'   (\code{custom.note = ""}), the note will be omitted completely.
#' @param override.coef Set custom values for the coefficients. New coefficients
#'   are provided as a list of numeric vectors. The list contains vectors of 
#'   coefficients for each model. There must be as many vectors of coefficients 
#'   as there are models. For example, if there are two models with three model 
#'   terms each, the argument could be specified as 
#'   \code{override.coef = list(c(0.1, 0.2, 0.3), c(0.05, 0.06, 0.07))}. 
#'   If there is only one model, custom values can be provided as a plain vector
#'   (not embedded in a list). 
#'   For example: \code{override.coef = c(0.05, 0.06, 0.07)}.
#' @param override.se Set custom values for the standard errors. This only 
#'   has an effect where standard errors are converted into confidence 
#'   intervals because no other CIs are present. New standard errors are 
#'   provided as a list of numeric vectors. The list contains vectors of 
#'   standard errors for each model. There must be as many vectors of 
#'   standard errors as there are models. For example, if there are two
#'   models with three coefficients each, the argument could be specified 
#'   as \code{override.se = list(c(0.1, 0.2, 0.3), c(0.05, 0.06, 0.07))}. 
#'   If there is only one model, custom values can be provided as a plain
#'   vector (not embedded in a list). For example: 
#'   \code{override.se = c(0.05, 0.06, 0.07)}. Overriding standard 
#'   errors can be useful for the implementation of robust SEs, for example.
#' @param override.pval Set custom values for the p values. This only has
#'   an effect where standard errors are converted into confidence intervals
#'   because no other CIs are present. In this case, significance is derived 
#'   from the p values rather than the confidence intervals. New p values are 
#'   provided as a list of numeric vectors. The list contains vectors of 
#'   p values for each model. There must be as many vectors of p values as 
#'   there are models. For example, if there are two models with three 
#'   coefficients each, the argument could be specified as 
#'   \code{override.pval = list(c(0.1, 0.2, 0.3), c(0.05, 0.06, 0.07))}. 
#'   If there is only one model, custom values can be provided as a plain 
#'   vector (not embedded in a list). For example: 
#'   \code{override.pval = c(0.05, 0.06, 0.07)}. 
#'   Overriding p values can be useful for the implementation 
#'   of robust SEs and p values, for example.
#' @param override.ci.low Set custom lower confidence interval bounds. 
#'   This works like the other override arguments, with one exception: if 
#'   confidence intervals are provided here and in the \code{override.ci.up} 
#'   argument, the standard errors and p values as well as the \code{ci.force} 
#'   argument are ignored.
#' @param override.ci.up Set custom upper confidence interval bounds. 
#'   This works like the other override arguments, with one exception: if 
#'   confidence intervals are provided here and in the \code{override.ci.low} 
#'   argument, the standard errors and p values as well as the \code{ci.force} 
#'   argument are ignored.
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
#'   function or \code{\link{extract-methods}}.
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
#' @param ci.level If standard errors are converted to confidence intervals
#'   (because a model does not natively support CIs), what confidence level 
#'   should be used for the outer confidence interval? By default, \code{0.95}
#'   is used (i.e., an alpha value of 0.05).
#' @param use.se Use one standard error for the inner horizontal bar and two 
#'   standard errors from the estimate for the outer horizontal bar (instead 
#'   of confidence intervals). Only available if standard errors can be 
#'   extracted from the model using the respective \code{\link{extract}} function.
#' @param type The default option is \code{\link{facet}}. If one model only is 
#'   specified, it will print one forest plot applied to point estimates and 
#'   confidence intervals. If more than one model is specified, it will print 
#'   as many plots as the number of models in a column of plots. Alternatively,
#'   if \code{\link{forest}} is specified, coefficients from one or more models 
#'   will be grouped together and displayed as a forest plot.
#' @param theme The \code{\link{theme}} argument can be used to customomise 
#'   the non data component of your plot. The default \code{\link{theme}} is
#'   \code{theme_bw} and it can be substituted by any other theme compatible 
#'   with \pkg{ggplot2}.
#' @param custom.coef.map The \code{custom.coef.map} argument can be used to
#'   select, omit, rename, and reorder coefficients.
#'   Users must supply a named list of this form: \code{list("x" = "First
#'   variable", "y" = NA, "z" = "Third variable")}. With that particular example
#'   of \code{custom.coef.map},
#'   \enumerate{
#'   \item coefficients will presented in order: \code{"x"}, \code{"y"},
#'   \code{"z"}.
#'   \item variable \code{"x"} will appear as \code{"First variable"}, variable
#'   \code{"y"} will appear as \code{"y"}, and variable \code{"z"} will
#'   appear as "Third variable".
#'   \item all variables not named \code{"x"}, \code{"y"}, or \code{"z"} will
#'   be omitted from the table.}
#' @param signif.light Color of outer confidence intervals for significant model terms.
#' @param signif.medium Color of inner confidence intervals for significant model terms.
#' @param signif.dark Color of point estimates and labels for significant model terms.
#' @param insignif.light Color of outer confidence intervals for insignificant model terms.
#' @param insignif.medium Color of inner confidence intervals for insignificant model terms.
#' @param insignif.dark Color of point estimates and labels for insignificant model terms.
#' @param ggsave.output Using this argument, the resulting table is written to a file 
#'   rather than to the R prompt. The file name can be specified as a character string. 
#'   The file extension is automatically recognized. \code{pdf}, \code{ps}, \code{png}, 
#'   \code{bmp}, \code{jpg}, and \code{tiff} are supported.
#' @param ... Custom options to be passed on to the extract function or the 
#'   gggplot2 device. See the help entries of \link{extract} and \link{extract-methods} 
#'   for more information. 
#' @return coefficient plots from statistical model output using \code{\link{ggplot2}}
#'
#' @author Philip Leifeld
#' @family texreg
#' @seealso \code{\link{texreg-package}} \code{\link{extract}}
#'   \code{\link{extract-methods}} \code{\link{texreg}}
#'   
#'   @examples
#' # example from the 'lm' help file:
#' ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#' trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#' group <- gl(2,10,20, labels = c("Ctl","Trt"))
#' weight <- c(ctl, trt)
#' lm.D9 <- lm(weight ~ group)
#' screenreg(lm.D9)  # print model output to the R console
#' plotreg(lm.D9)    # plot model output as a diagram
#' # plot model output as a diagram customising theme and title and automatically saving in pdf.
#' dev.off()
#' plotreg(lm.D9, theme = theme_dark(),  ggtitle ="my title", ggsave.output = "myplot.pdf")
#' unlink("myplot.pdf")
#' 
#' @export
plotreg <- function(l,
                    custom.model.names = NULL,
                    custom.title = NULL,
                    custom.coef.names = NULL,
                    custom.note = NULL,
                    override.coef = 0,
                    override.se = 0,
                    override.pval = 0,
                    override.ci.low = 0,
                    override.ci.up = 0,
                    omit.coef = NULL,
                    reorder.coef = NULL,
                    ci.level = 0.95,
                    use.se = FALSE,
                    type = "facet",
                    theme = theme_bw(),
                    custom.coef.map = NULL,
                    signif.light = "#FBC9B9",
                    signif.medium = "#F7523A",
                    signif.dark = "#BD0017",
                    insignif.light = "#C5DBE9",
                    insignif.medium = "#5A9ECC",
                    insignif.dark = "#1C5BA6",
                    ggsave.output = NULL,
                    ...) {
    
    if (!is.null(omit.coef) && !is.na(omit.coef) && !is.character(omit.coef)) {
        stop("omit.coef must be a character string!")
    }
    if (is.null(omit.coef)) {
        omit.coef <- NA
    }
    
    if (length(use.se) == 1) {
        use.se <- rep(use.se, length(l))
    }
    
    # extract texreg objects and override data
    models <- get.data(l, ...)
    models <- override(models, override.coef, override.se, override.pval, 
                       override.ci.low, override.ci.up)
    # custom model names
    model.names <- character()
    if (is.null(custom.model.names)) {
        model.names <- paste("Model", 1:length(l))
    } else if (length(custom.model.names) == 1) {
        model.names <- rep(custom.model.names, length(l))
    } else if (length(custom.model.names) != length(l)) {
        stop(paste("The 'custom.model.names' argument must have the same length as", "the 'l' argument."))
    } else {
        model.names <- custom.model.names
    }
    
    co <- se <- co.names <- pv <- lab <- ci.low <- ci.up <- NULL
    
    # make dataframe
    for (i in 1:length(models)) {
        # custom coefficient names
        if (!is.null(custom.coef.names)) {
            if (class(custom.coef.names) == "list") {
                if (length(custom.coef.names[[i]]) == length(models[[i]]@coef.names)) {
                    models[[i]]@coef.names <- custom.coef.names[[i]]
                } else {
                    stop(paste0("Model ", i, ": wrong number of custom coefficient names."))
                }
            } else if (class(custom.coef.names) == "character") {
                if (length(custom.coef.names) == length(models[[i]]@coef.names)) {
                    models[[i]]@coef.names <- custom.coef.names
                } else {
                    stop(paste0("Model ", i, ": wrong number of custom coefficient names."))
                }
            } else {
                stop(paste("Custom coefficient names must be provided as a list of", "character vectors."))
            }
        }
        
        coef.names <- models[[i]]@coef.names
        co.names <- append(co.names, coef.names)
        coefs <- models[[i]]@coef
        co <- append(co, coefs)
        label <- rep(paste0("Model", i), length(models[[i]]@coef))
        lab <- append(lab, label)
        if (length(models[[i]]@se) > 0) {
            s.e <- models[[i]]@se
            se <- append(se, s.e)
        } 
        if (length(models[[i]]@pvalues) > 0) {
            pvalues <- models[[i]]@pvalues
            pv <- append(pv, pvalues)
        }
        if (length(models[[i]]@ci.low) > 0) {
            ci.lower <- models[[i]]@ci.low
            ci.low <- append(ci.low, ci.lower)
        }
        if (length(models[[i]]@ci.up) > 0) {
            ci.upper <- models[[i]]@ci.up
            ci.up <- append(ci.up, ci.upper)
        }
    }    
    dataframe <- data.frame(cbind(co.names, co, lab))
    
    if (length(models[[i]]@se) > 0) {
        dataframe <- cbind(dataframe, se)
    } 
    
    if (length(models[[i]]@pvalues) > 0) {
        dataframe <- cbind(dataframe, pv)
    } 
    
    if (length(models[[i]]@ci.low) > 0) {
        dataframe <- cbind(dataframe, ci.low)
    } 
    
    if (length(models[[i]]@ci.up) > 0) {
        dataframe <- cbind(dataframe, ci.up)
    } 
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
    
    # aggregate confidence intervals or SEs
    if (isTRUE(use.se[i]) && length(models[[i]]@se) > 0) {
        lower.inner <- dataframe$co - dataframe$se
        upper.inner <- dataframe$co + dataframe$se
        lower.outer <- dataframe$co - 2 * dataframe$se
        upper.outer <- dataframe$co + 2 * dataframe$se
    } else if (length(models[[i]]@ci.low) == 0 && length(models[[i]]@se) > 0) {
        z.inner <- qnorm(1 - ((1 - 0.5) / 2))
        lower.inner <- dataframe$co - (z.inner * dataframe$se)
        upper.inner <- dataframe$co + (z.inner * dataframe$se)
        z.outer <- qnorm(1 - ((1 - ci.level) / 2))
        lower.outer <- dataframe$co - (z.outer * dataframe$se)
        upper.outer <- dataframe$co + (z.outer * dataframe$se)
    } else if (length(models[[i]]@se) == 0 && length(models[[i]]@pvalues) > 0) {
        stop("Model has p-values but no SEs. SEs or CIs are required for plotting.")
    } else {
        lower.outer <- dataframe$ci.low
        upper.outer <- dataframe$ci.up
        lower.inner <- rep(0, length(dataframe$co.names))
        upper.inner <- rep(0, length(dataframe$co.names))
    }
    
    if (length(lower.outer) > 0) {
        dataframe <- cbind(dataframe, lower.outer)     
    }
    if (length(lower.inner) > 0) {
        dataframe <- cbind(dataframe, lower.inner)     
    }
    if (length(upper.outer) > 0) {
        dataframe <- cbind(dataframe, upper.outer)     
    }    
    if (length(upper.inner) > 0) {
        dataframe <- cbind(dataframe, upper.inner)     
    }
    
    # # which terms are significant?
    if (length(dataframe$pv) > 0) {
        signif.outer <- dataframe$pv < (1 - ci.level)
    } else if (length(dataframe$ci.low) > 0) {
        signif.outer <- ((dataframe$ci.low) > 0 & (dataframe$ci.up) > 0 | (dataframe$ci.low) < 0 & (dataframe$ci.up) < 0)
    }
    
    if (sapply(signif.outer, is.logical) &&
        length(signif.outer) == length(dataframe$co)) {
        signif <- sapply(signif.outer, isTRUE)
    } else if (sapply(signif.outer, isTRUE)) {
        signif <- apply(cbind(lower.outer, upper.outer), 1,
                        function(x) x[1] <= 0 && x[2] >= 0)
    } else if (sapply(signif.outer, isFALSE)) {
        signif <- apply(cbind(lower.inner, upper.inner), 1,
                        function(x) x[1] <= 0 && x[2] >= 0)
    } else {
        stop("signif.outer does not correspond to the intervals provided.")
    }
    
    
    dataframe <- cbind(dataframe, signif)
    
    if (length(co) == 0) {
        stop(paste("No coefficients available. Was the 'omit.coef' argument", 
                   "misspecified? If coefficients were renamed using the",
                   "'custom.coef.names' argument, 'omit.coef' refers to the renamed",
                   "coefficients."))
    }
    
    # omit.coef
    if (!is.na(omit.coef)) {
        dataframe <- dataframe[!grepl(omit.coef, dataframe$co.names), ]
    }
    
    # custom coef. name
    if (length(custom.coef.names) > 0) {
        levels(dataframe$co.names) <- custom.coef.names
    }
    
    # custom model name
    if (length(custom.model.names) > 0) {
        levels(dataframe$lab) <- custom.model.names
    }
    
    
    # ggplot functions   
    if (type == "facet") {
        
        if (length(reorder.coef) > 0) {
            # reorder levels 
            dataframe$co.names <- ordered(dataframe$co.names, levels = (unique(as.character(dataframe$co.names))))
            # reorder coef
            reorder.coef<- rev(reorder.coef)
            dataframe$co.names <- factor(dataframe$co.names, levels(dataframe$co.names)[reorder.coef])
        } else {
            # reorder levels 
            dataframe$co.names <- ordered(dataframe$co.names, levels = (unique(as.character(dataframe$co.names))))
            # reorder coef in original order
            startlevel<- c(length(dataframe$co.names):1)
            dataframe$co.names <- factor(dataframe$co.names, levels(dataframe$co.names)[startlevel])
        }
        
        if (length(custom.coef.map) > 0) {
            # keep only coefficients selected by the user and replace NAs if any
            idx <- is.na(custom.coef.map)
            custom.coef.map[idx] <- names(custom.coef.map)[idx]
            selectcoef <- names(custom.coef.map)
            keepcoef <- paste(selectcoef, collapse="|")
            dataframe <- dataframe[grepl(keepcoef, dataframe$co.names), ]
            # resize number of levels in factor to avoid error
            dataframe$co.names <- ordered(dataframe$co.names, levels = (unique(as.character(dataframe$co.names))))
            # rename selected coefficients
            renamedcoef <- paste(unlist(custom.coef.map), sep = ", ")
            levels(dataframe$co.names) <- renamedcoef
            # invert order factors for plot
            coeforder<- c(length(dataframe$co.names):1)
            dataframe$co.names <- factor(dataframe$co.names, levels(dataframe$co.names)[coeforder])
        }
        
        p <- ggplot(dataframe, aes(co.names, co)) + 
            geom_hline(yintercept = 0, 
                       lty = 2, lwd = 1, 
                       colour="grey50") +
            geom_errorbar(aes(ymin=lower.outer, ymax=upper.outer), 
                          lwd = 1, 
                          colour = ifelse(sapply(dataframe$signif, isTRUE) , signif.light, insignif.light), 
                          width = 0) +
            geom_errorbar(aes(ymin= lower.inner, ymax= upper.inner), 
                          lwd = 2.5, 
                          colour = ifelse(sapply(dataframe$signif, isTRUE), signif.medium, insignif.medium), 
                          width = 0) +
            geom_point(size = 3, 
                       pch = ifelse(sapply(dataframe$signif, isTRUE), 21, 22), 
                       fill = ifelse(sapply(dataframe$signif, isTRUE), signif.dark, insignif.dark)) +
            coord_flip() 
        p <- p + theme
        p <- p + xlab(" ") +
            facet_wrap(~lab, 
                       strip.position="left", 
                       nrow=length(dataframe), 
                       scales = "free_y")
        
        if (length(custom.title) > 0) {
            p <- p + ggtitle(custom.title)
        } 
        
    } else if (type == "forest") {
        if (length(reorder.coef) > 0) {
            # reorder levels 
            dataframe$co.names <- ordered(dataframe$co.names, levels = (unique(as.character(dataframe$co.names))))
            # reorder coef
            dataframe$co.names <- factor(dataframe$co.names, levels(dataframe$co.names)[reorder.coef])
        } else {
            # reorder levels 
            dataframe$co.names <- ordered(dataframe$co.names, levels = (unique(as.character(dataframe$co.names))))
            #reorder coef in original order
            startlevel<- c(1:length(dataframe$co.names))
            dataframe$co.names <- factor(dataframe$co.names, levels(dataframe$co.names)[startlevel])
        }
        
        if (length(custom.coef.map) > 0) {
            # keep only coefficients selected by the user and replace NAs if any
            idx <- is.na(custom.coef.map)
            custom.coef.map[idx] <- names(custom.coef.map)[idx]
            selectcoef <- names(custom.coef.map)
            keepcoef <- paste(selectcoef, collapse="|")
            dataframe <- dataframe[grepl(keepcoef, dataframe$co.names), ]
            # resize number of levels in factor to avoid error
            dataframe$co.names <- ordered(dataframe$co.names, levels = (unique(as.character(dataframe$co.names))))
            # rename selected coefficients
            renamedcoef <- paste(unlist(custom.coef.map), sep = ", ")
            levels(dataframe$co.names) <- renamedcoef
            # invert order factors for plot
            coeforder<- c(length(dataframe$lab):1)
            dataframe$lab <- factor(dataframe$lab, levels(dataframe$lab)[coeforder])
        }
        
        # put models in the right order
        laborder<- c(length(dataframe$lab):1)
        dataframe$lab <- factor(dataframe$lab, levels(dataframe$lab)[laborder])
        
        # plot
        p <- ggplot(dataframe, aes(lab, co)) + 
            geom_hline(yintercept = 0, 
                       lty = 2, 
                       lwd = 1, 
                       colour="grey50") +
            geom_errorbar(aes(ymin=lower.outer, ymax=upper.outer), 
                          lwd = 1, 
                          colour = ifelse(sapply(dataframe$signif, isTRUE) , signif.light, insignif.light), 
                          width = 0) +
            geom_errorbar(aes(ymin= lower.inner, ymax= upper.inner), 
                          lwd = 2.5, 
                          colour = ifelse(sapply(dataframe$signif, isTRUE), signif.medium, insignif.medium),
                          width = 0) +
            geom_point(size = 3, 
                       pch = ifelse(sapply(dataframe$signif, isTRUE), 21, 22), 
                       fill = ifelse(sapply(dataframe$signif, isTRUE), signif.dark, insignif.dark)) +
            coord_flip() 
        p <- p + theme
        p <- p + xlab(" ") +
            facet_wrap(~co.names, 
                       strip.position="left", 
                       nrow=length(dataframe), 
                       scales = "free_y")
        
        if (length(custom.title) > 0) {
            p <- p + ggtitle(custom.title)
        } 
    } 
    
    # adds message to p as ylab; adds output meassages to paste that comes as R output (not in plot)
    if (isTRUE(use.se[i]) && length(models[[i]]@se) > 0) {
        p <- p + ylab("Bars denote SEs. Circle points denote significance.")
        if (length(custom.note) > 0) {
            p <- p + ylab(custom.note)
        }
        if (length(models) == 1) {
            message("Model: bars denote one (inner) resp. two (outer) standard errors.")
            
        } else if (length(models) > 1) {
            message("Models: bars denote one (inner) resp. two (outer) standard errors.")
        }   
    } else if (length(models[[i]]@ci.low) == 0 && length(models[[i]]@se) > 0) {
        p <- p + ylab("Bars denote CIs. Circle points denote significance.")
        if (length(custom.note) > 0) {
            p <- p + ylab(custom.note)
        }
        if (length(models) == 1) {
            message(paste0("Model: bars denote 0.5 (inner) resp. ", ci.level, 
                           " (outer) confidence intervals (computed from standard errors)."))
        } else if (length(models) > 1) {
            message(paste0("Models: bars denote 0.5 (inner) resp. ", ci.level, 
                           " (outer) confidence intervals (computed from standard errors)."))
        }
    } else {
        p <- p + ylab("Bars denote CIs. Circle points denote significance.")
        if (length(custom.note) > 0) {
            p <- p + ylab(custom.note)
        }
        if (length(models) == 1) {
            message(paste0("Model: bars denote  ", ci.level, " confidence intervals."))
        } else if (length(models) > 1) {
            message(paste0("Models: bars denote  ", ci.level, " confidence intervals."))
        }  
    }
    
    # print plot in R
    print(p)
    
    # file output - save plot with ggsave
    if (length(ggsave.output) > 0 ) {
        ggsave(ggsave.output, plot = p) 
    }
}
