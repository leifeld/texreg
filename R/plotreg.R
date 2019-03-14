# plotreg function
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
            dataframe <- custommap(dataframe, custom.coef.map)
        }
        
        # plot
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
            p <- p +   ggtitle(custom.title)
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
            #cat("\ndataframe1: ", deparse(dataframe))
            #dataframe <- custommap(dataframe, custom.coef.map)
            #cat("\ndataframe2: ", rownames(dataframe))
            
            # separate data frame by model
            rownames(dataframe) <- dataframe$co.names
            dataframes <- split(dataframe, dataframe$lab)
            
            # apply custommap to each model
            rownames(dataframe) <- dataframe$co.names
            lapply(dataframes, function x rownames())
            dataframes <- lapply(dataframes, custommap)
            cat("\ndataframes", deparse(dataframes))
            
            # putting them together again
            dataframes <- Reduce(rbind, dataframes)
            
            # when user supplies NA as destination, replace with origin
           # idx <- is.na(custom.coef.map)
           # custom.coef.map[idx] <- names(custom.coef.map)[idx]
            
            # subset of coefficients to keep
            # origins <- NULL
            # for (i in 1:length(dataframes)) {
            #     origin <- names(custom.coef.map)[names(custom.coef.map) %in% dataframe$co.names[i]]
            #     origins <- append(origins, rbind(origin))
            # }
            # cat("\norigins: ", origins)
            # destination <- unlist(custom.coef.map[origins])
            # # otherwise R converts to numeric if a single coefficient is passed
            # out <- m[origin, , drop = FALSE]
            # 
            # # rename
            # row.names(out) <- destination
            # 
            # # output
            # return(out)
            
            
            #dataframes <- lapply(dataframes, function(x) custommap(x, custom.coef.map))
            #cat("\nlength(dataframes): ", length(dataframes))
            #cat("\ndataframe1: ", deparse(dataframes))
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
            p <- p +   ggtitle(custom.title)
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
    
    
    if (!is.null(file) && !is.na(file)) {
        dev.off()
    }
}

