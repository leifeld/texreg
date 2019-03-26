#' Create coefficient plots from statistical model output using \pkg{ggplot2}.
#'
#' Create coefficient plots of \R regression output using \pkg{ggplot2}.
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
#'  possible to display the plots in two ways: using the \code{type = "facet"} 
#'  option, one forest plot applied to point estimates and confidence intervals 
#'  will be visualised in case there is only one model. If there is more than 
#'  one model each one will be plotted next to the other; using the 
#'  \code{type = "forest"} option, coefficients from one or more models will 
#'  be grouped together and displayed as a forest plot. 
#'
#' @param override.pval Replacement list of p-value vectors.
#' @param ci.level If standard errors are converted to confidence intervals
#'   (because a model does not natively support CIs), what confidence level 
#'   should be used for the outer confidence interval? By default, \code{0.95}
#'   is used (i.e., an alpha value of 0.05).
#' @param type The default option is \code{type = "facet"}. If one model only is 
#'   specified, it will print one forest plot applied to point estimates and 
#'   confidence intervals. If more than one model is specified, it will print 
#'   as many plots as the number of models in a column of plots. Alternatively,
#'   if \code{type = "forest"} is specified, coefficients from one or more models 
#'   will be grouped together and displayed as a forest plot.
#' @param theme The \code{\link{theme}} argument can be used to customomise 
#'   the non data component of your plot. The default \code{\link{theme}} is
#'   \code{theme_bw} and it can be substituted by any other theme compatible 
#'   with \pkg{ggplot2}.
#' @param signif.light Color of outer confidence intervals for significant model terms.
#' @param signif.medium Color of inner confidence intervals for significant model terms.
#' @param signif.dark Color of point estimates and labels for significant model terms.
#' @param insignif.light Color of outer confidence intervals for insignificant model terms.
#' @param insignif.medium Color of inner confidence intervals for insignificant model terms.
#' @param insignif.dark Color of point estimates and labels for insignificant model terms.
#' 
#' @return coefficient plots from statistical model output using \pkg{ggplot2}.
#' 
#' @inheritParams texreg
#' @inheritParams matrixreg
#'
#'
#' @author Philip Leifeld
#' @family texreg
#' @seealso \code{\link{texreg-package}} \code{\link{extract}}
#'   \code{\link{extract-methods}} \code{\link{texreg}}
#'   
#' @examples
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
                    file = NULL,
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
                    ...) {
    
    if (!is.null(omit.coef) && !is.na(omit.coef) && !is.character(omit.coef)) {
        stop("'omit.coef' must be a character string!")
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
    
    # which terms are significant?
    if (length(dataframe$pv) > 0) {
        signif.outer <- dataframe$pv < (1 - ci.level)
    } else if (length(dataframe$ci.low) > 0) {
        signif.outer <- ((dataframe$ci.low) > 0 & (dataframe$ci.up) > 0 | (dataframe$ci.low) < 0 & (dataframe$ci.up) < 0)
    }
    
    if (is.logical(signif.outer) && length(signif.outer) == length(dataframe$co)) {
        signif <- sapply(signif.outer, isTRUE)
    } else if (sapply(signif.outer, isTRUE)) {
        signif <- apply(cbind(lower.outer, upper.outer), 1,
                        function(x) x[1] <= 0 && x[2] >= 0)
    } else if (sapply(signif.outer, isFALSE)) {
        signif <- apply(cbind(lower.inner, upper.inner), 1,
                        function(x) x[1] <= 0 && x[2] >= 0)
    } else {
        stop("'signif.outer' does not correspond to the intervals provided.")
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
            keepcoef <- paste(selectcoef, collapse = "|")
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
                       colour = "grey50") +
            geom_errorbar(aes(ymin = lower.outer, ymax = upper.outer), 
                          lwd = 1, 
                          colour = ifelse(sapply(dataframe$signif, isTRUE) , signif.light, insignif.light), 
                          width = 0) +
            geom_errorbar(aes(ymin = lower.inner, ymax = upper.inner), 
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
                       strip.position = "left", 
                       nrow = length(dataframe), 
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
            # reorder coef in original order
            startlevel<- c(1:length(dataframe$co.names))
            dataframe$co.names <- factor(dataframe$co.names, levels(dataframe$co.names)[startlevel])
        }
        
        if (length(custom.coef.map) > 0) {
            # keep only coefficients selected by the user and replace NAs if any
            idx <- is.na(custom.coef.map)
            custom.coef.map[idx] <- names(custom.coef.map)[idx]
            selectcoef <- names(custom.coef.map)
            keepcoef <- paste(selectcoef, collapse = "|")
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
                       colour = "grey50") +
            geom_errorbar(aes(ymin = lower.outer, ymax = upper.outer), 
                          lwd = 1, 
                          colour = ifelse(sapply(dataframe$signif, isTRUE) , signif.light, insignif.light), 
                          width = 0) +
            geom_errorbar(aes(ymin = lower.inner, ymax = upper.inner), 
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
                       strip.position = "left", 
                       nrow = length(dataframe), 
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
    if (is.null(file)) {
        return(p)
    } else {
        ggsave(file = file, plot = p)
    }
}
