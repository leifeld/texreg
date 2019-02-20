# # function which plots coefficients and confidence intervals ('coefplot')
# #coefplot <- function(labels, estimates, lower.inner = NULL, 
#     upper.inner = NULL, lower.outer = NULL, upper.outer = NULL, 
#     signif.outer = TRUE, xlab = "Coefficients and confidence intervals", 
#     main = "Coefficient plot", xlim = NULL, cex = 2.5, lwd.zerobar = 4, 
#     lwd.vbars = 1, lwd.inner = 7, lwd.outer = 5, ylab.cex = 1.0, 
#     signif.light = "#fbc9b9", signif.medium = "#f7523a", 
#     signif.dark = "#bd0017", insignif.light = "#c5dbe9", 
#     insignif.medium = "#5a9ecc", insignif.dark = "#1c5ba6", ...) {
#   
#   # check consistency of arguments
#   if (length(table(c(length(estimates), length(lower.outer), 
#       length(upper.outer), length(labels)))) > 1) {
#     stop("Vectors should have the same length.")
#   }
#   
#   # define shortcut variables outer, mn, mx, steps, and num
#   if (!is.null(lower.outer) && !is.null(upper.outer)) {
#     outer <- TRUE
#   } else {
#     outer <- FALSE
#   }
#   if (!is.null(lower.inner) && !is.null(upper.inner)) {
#     inner <- TRUE
#   } else {
#     inner <- FALSE
#   }
#   if (outer == TRUE) {
#     mn <- min(lower.outer)
#     mx <- max(upper.outer)
#   } else if (inner == TRUE) {
#     mn <- min(lower.inner)
#     mx <- max(upper.inner)
#   } else {
#     stop("Either inner or outer bounds must be provided.")
#   }
#   if (mn > 0) {
#     mn <- 0
#   }
#   if (mx < 0) {
#     mx <- 0
#   }
#   if (!is.null(xlim) && length(xlim) == 2 && is.numeric(xlim)) {
#     mn <- xlim[1]
#     mx <- xlim[2]
#   }
#   steps <- floor(mn):ceiling(mx)
#   num <- length(labels)
#   
#   
#   # which terms are significant?
#   if (class(signif.outer) == "logical" && 
#       length(signif.outer) == length(estimates)) {
#     signif <- signif.outer == FALSE
#   } else if (outer == TRUE && signif.outer == TRUE) {
#     signif <- apply(cbind(lower.outer, upper.outer), 1, 
#         function(x) x[1] <= 0 && x[2] >= 0)
#   } else if (inner == TRUE && signif.outer == FALSE) {
#     signif <- apply(cbind(lower.inner, upper.inner), 1, 
#         function(x) x[1] <= 0 && x[2] >= 0)
#   } else {
#     stop("signif.outer does not correspond to the intervals provided.")
#   }
#   signif <- signif == FALSE
#   
#   # create plot; compute left margin; add labels to left margin
#   inches.to.lines <- (par("mar") / par("mai"))[1]
#   lab.width <- max(strwidth(labels, units = "inches") * inches.to.lines * ylab.cex)
#   
#   par(mar = c(5, 1 + lab.width, 4, 2) + 0.1)
#   plot(0, 0, xlim = c(mn, mx), ylim = c(1, length(labels)), xlab = xlab, 
#       main = main, ylab = "", axes = FALSE, type = "n", ...)
#   axis(side = 2, at = 1:num, labels = FALSE, las = 2, 
#       tck = 0, lty = 0, ...)
#   if (any(signif)) {
#     mtext(side = 2, text = labels[which(signif)], at = num + 1 - which(signif), 
#         col = signif.dark, las = 2, cex = ylab.cex, ...)
#   }
#   if (any(!signif)) {
#     mtext(side = 2, text = labels[which(!signif)], 
#         at = num + 1 - which(!signif), col = insignif.dark, las = 2, 
#         cex = ylab.cex, ...)
#   }
#   axis(side = 1, las = 0, lty = 0, ...)
#   
#   # add vertical gray lines
#   zeros <- rep(0, length(steps))
#   ends <- rep((num + 1), length(steps))
#   segments(steps, zeros, steps, ends, col = "gray", lwd = lwd.vbars)
#   segments(0, 0, 0, num + 1, lwd = lwd.zerobar, col = "grey75")
#   
#   # draw outer CIs
#   if (outer == TRUE) {
#     if (any(signif)) {
#       segments(lower.outer[signif], num + 1 - which(signif), 
#           upper.outer[signif], num + 1 - which(signif), lwd = lwd.outer, 
#           col = signif.light)
#     }
#     if (any(!signif)) {
#       segments(lower.outer[!signif], num + 1 - which(!signif), 
#           upper.outer[!signif], num + 1 - which(!signif), lwd = lwd.outer, 
#           col = insignif.light)
#     }
#   }
#   
#   # draw inner CIs
#   if (inner == TRUE) {
#     if (any(signif)) {
#       segments(lower.inner[signif], num + 1 - which(signif), 
#           upper.inner[signif], num + 1 - which(signif), lwd = lwd.inner, 
#           col = signif.medium)
#     }
#     if (any(!signif)) {
#       segments(lower.inner[!signif], num + 1 - which(!signif), 
#           upper.inner[!signif], num + 1 - which(!signif), lwd = lwd.inner, 
#           col = insignif.medium)
#     }
#     if (outer == FALSE) {
# 
#     }
#   }
#   
#   # draw point estimates
#   if (any(signif)) {
#     points(estimates[signif], num + 1 - which(signif), col = signif.dark, 
#         pch = 20, cex = cex)
#   }
#   if (any(!signif)) {
#     points(estimates[!signif], num + 1 - which(!signif), col = insignif.dark, 
#         pch = 20, cex = cex)
#   }
# }
# 

# plotreg function
plotreg <- function(l, file = NULL, custom.model.names = NULL, 
                    custom.coef.names = NULL, custom.note = NULL, override.coef = 0, 
                    override.se = 0, override.pval = 0, override.ci.low = 0, 
                    override.ci.up = 0, omit.coef = NULL, reorder.coef = NULL, 
                    ci.level = 0.95, use.se = FALSE, mfrow = TRUE, xlim = NULL, 
                    cex = 2.5, lwd.zerobar = 4, lwd.vbars = 1, lwd.inner = 7, 
                    lwd.outer = 5, ylab.cex = 1.0, signif.light = "#fbc9b9", 
                    signif.medium = "#f7523a", signif.dark = "#bd0017", 
                    insignif.light = "#c5dbe9", insignif.medium = "#5a9ecc", 
                    insignif.dark = "#1c5ba6", type = "facet", co = NULL, se = NULL,
                    co.names = NULL, pv = NULL, lab = NULL, ci.low = NULL, ci.up = NULL, ...) {
    
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
        stop(paste("The 'custom.model.names' argument must have the same length as",
                   "the 'l' argument."))
    } else {
        model.names <- custom.model.names
    }
    
    # file output; check file extension
    if (!is.null(file) && !is.na(file)) {
        if (grepl(".pdf$", file)) {
            pdf(file, ...)
        } else if (grepl(".jpe*g$", file)) {
            jpeg(file, ...)
        } else if (grepl(".png$", file)) {
            png(file, ...)
        } else if (grepl(".bmp$", file)) {
            bmp(file, ...)
        } else if (grepl(".tiff*$", file)) {
            tiff(file, ...)
        } else if (grepl(".ps$", file)) {
            postscript(file, ...)
        } else {
            stop("File extension not recognized.")
        }
    }
    
    # mfrow argument: arrange multiple plots on a page
    if (mfrow == TRUE) {
        if (length(models) == 1) {
            par(mfrow = c(1, 1))
        } else if (length(models) == 2) {
            par(mfrow = c(2, 1))
        } else if (length(models) == 3 || length(models) == 4 || 
                   length(models) > 6) {
            par(mfrow = c(2, 2))
        } else if (length(models) == 5 || length(models) == 6) {
            par(mfrow = c(3, 2))
        }
    }
    
    # plot each model
    for (i in 1:length(models)) {
        
        # custom coefficient names
        if (!is.null(custom.coef.names)) {
            if (class(custom.coef.names) == "list") {
                if (length(custom.coef.names[[i]]) == length(models[[i]]@coef.names)) {
                    models[[i]]@coef.names <- custom.coef.names[[i]]
                } else {
                    stop(paste0("Model ", i, 
                                ": wrong number of custom coefficient names."))
                }
            } else if (class(custom.coef.names) == "character") {
                if (length(custom.coef.names) == length(models[[i]]@coef.names)) {
                    models[[i]]@coef.names <- custom.coef.names
                } else {
                    stop(paste0("Model ", i, 
                                ": wrong number of custom coefficient names."))
                }
            } else {
                stop(paste("Custom coefficient names must be provided as a list of", 
                           "character vectors."))
            }
        }
        
        # remove any of the coefficients? apply regex
        #remove.rows <- grep(omit.coef, models[[i]]@coef.names)
        #if (is.na(omit.coef)) {
        #    remove.rows <- length(models[[i]]@coef) + 1
        #}
        
        coef.names <- models[[i]]@coef.names
        co.names <- append(co.names, coef.names)
        coef <- models[[i]]@coef
        co <- append(co, coef)
        label <- rep(paste0("Model", i),length(models[[i]]@coef))
        lab <- append(lab, label)
        if (length(models[[i]]@se) > 0) {
            s.e <- models[[i]]@se
            se <- append(se, s.e)
        } 
        if(length(models[[i]]@pvalues) > 0) {
            pvalues <- models[[i]]@pvalues
            pv <- append(pv, pvalues)
        }
        if(length(models[[i]]@ci.low) > 0) {
            ci.lower <- models[[i]]@ci.low
            ci.low <- append(ci.low, ci.lower)
        }
        if(length(models[[i]]@ci.up) > 0) {
            ci.upper <- models[[i]]@ci.up
            ci.up <- append(ci.up, ci.upper)
        }
        
    }    
    dataframe <- data.frame(cbind(co.names, co, lab ))
    
    if (length(models[[i]]@se) > 0) {
        dataframe <- cbind(dataframe, se)
    } else {
        #do  nothing
    }
    if (length(models[[i]]@pvalues) > 0) {
        dataframe <- cbind(dataframe, pv)
    } else {
        #do  nothing
    }
    if (length(models[[i]]@ci.low) > 0) {
        dataframe <- cbind(dataframe, ci.low)
    } else {
        #do  nothing
    }
    if (length(models[[i]]@ci.up) > 0) {
        dataframe <- cbind(dataframe, ci.up)
    } else {
        #do  nothing
    }
    
    #reorder 
    dataframe <- reorder(dataframe, reorder.coef)
    
    #omit.coef
    #if (length(omit.coef)> 0){
    #    dataframe <- dataframe[grep(omit.coef, dataframe$co.names, invert = TRUE), ]
    #}
    #cat("\ndataframe: ", deparse(dataframe))
    dataframe$co <- as.numeric(as.character(dataframe$co))
    
    if (length(dataframe$se) > 0) {
        dataframe$se <- as.numeric(as.character(dataframe$se))
    } else {
        #do  nothing
    }
    if (length(dataframe$pv) > 0) {
        dataframe$pv <- as.numeric(as.character(dataframe$pv))
    } else {
        #do  nothing
    }
    if (length(dataframe$ci.low) > 0) {
        dataframe$ci.low <- as.numeric(as.character(dataframe$ci.low))
    } else {
        #do  nothing
    }
    if (length(dataframe$ci.up) > 0) {
        dataframe$ci.up <- as.numeric(as.character(dataframe$ci.up))
    } else {
        #do  nothing
    }
    cat("\ndataframe: ", deparse(dataframe))
    
    # aggregate confidence intervals or SEs
    if (use.se[i] == TRUE && length(models[[i]]@se) > 0) {
        lower.inner <- dataframe$co - dataframe$se
        upper.inner <- dataframe$co + dataframe$se
        lower.outer <- dataframe$co - 2 * dataframe$se
        upper.outer <- dataframe$co + 2 * dataframe$se
        # if (length(models[[i]]@pvalues) > 0) {
        #     dataframe$pv <- dataframe$pv[-remove.rows]
        #     dataframe$signif.outer <- dataframe$pv < (1 - ci.level)
        
        # sort according to reorder.coef argument
        #     dataframe <- data.frame(co.names, co, se, pv)
        #     dataframe <- reorder(dataframe, reorder.coef)
        #     pv <- dataframe[, 4]
        #} else {
        # sort according to reorder.coef argument
        #    dataframe <- data.frame(co.names, co, se)
        #    dataframe <- reorder(dataframe, reorder.coef)
        #}
        
    } else if (length(models[[i]]@ci.low) == 0 && length(models[[i]]@se) > 0) {
        z.inner <- qnorm(1 - ((1 - 0.5) / 2))
        lower.inner <- dataframe$co - (z.inner * dataframe$se)
        upper.inner <- dataframe$co + (z.inner * dataframe$se)
        z.outer <- qnorm(1 - ((1 - ci.level) / 2))
        lower.outer <- dataframe$co - (z.outer * dataframe$se)
        upper.outer <- dataframe$co + (z.outer * dataframe$se)
        #signif.outer <- TRUE
    } else if (length(models[[i]]@se) == 0 && length(models[[i]]@pvalues) > 0) {
        stop("Model has p-values but no SEs. SEs or CIs are required for plotting.")
    } else {
        #si può eliminare assegnando al data frame sopra direttamente lower.outer upper.outer
        lower.outer <- dataframe$ci.low
        upper.outer <- dataframe$ci.up
        lower.inner <- rep(0, length(dataframe))
        upper.inner <- rep(0, length(dataframe))
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
    
    if (length(dataframe$pv)> 0) {
        signif.outer <- dataframe$pv < (1 - ci.level)
    } 
    
    if(length(dataframe$ci.low) >0) {
        signif.outer <- ((dataframe$ci.low > 0 & dataframe$ci.up) > 0 | (dataframe$ci.low < 0 & dataframe$ci.up) < 0)
    }
    
    # which terms are significant?
    if (class(signif.outer) == "logical" &&
        length(signif.outer) == length(dataframe$co)) {
        signif <- signif.outer == FALSE
    } else if (outer == TRUE && signif.outer == TRUE) {
        signif <- apply(cbind(lower.outer, upper.outer), 1,
                        function(x) x[1] <= 0 && x[2] >= 0)
    } else if (inner == TRUE && signif.outer == FALSE) {
        signif <- apply(cbind(lower.inner, upper.inner), 1,
                        function(x) x[1] <= 0 && x[2] >= 0)
    } else {
        stop("signif.outer does not correspond to the intervals provided.")
    }
    signif <- signif == FALSE
    
    dataframe <- cbind(dataframe, signif)
    
    #place ggplot function below  
    if (type == "facet") {
        p <- ggplot(dataframe, aes(co.names, co)) + 
            geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
            geom_errorbar(aes(ymin=lower.outer, ymax=upper.outer), 
                          lwd=1, colour = ifelse(signif == TRUE , signif.light, insignif.light), width=0) +
            geom_errorbar(aes(ymin= lower.inner, ymax= upper.inner), lwd=2.5, 
                          colour = ifelse(signif == TRUE, signif.medium, insignif.medium), width=0) +
            geom_point(size=3, pch = ifelse(signif == TRUE, 21, 22), fill = ifelse(signif == TRUE, signif.dark, insignif.dark)) +
            coord_flip() + 
            theme_bw()
        if (length(models) == 1) {
            #cat(length(models))
            p <- p +   ggtitle("Model")
        } else if (length(models) > 1) {
            p <- p +  facet_wrap(~lab, strip.position="left", nrow=length(dataframe), scales = "free_y") +
                ggtitle("Models")
        }
        
    } else if (type == "forest") {
        cat("\ndataframeforest: ", deparse(dataframe))
        p <- ggplot(dataframe, aes(lab, co)) + 
            geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
            geom_errorbar(aes(ymin= lower.inner, ymax= upper.inner), 
                          lwd=1, colour = ifelse(signif == TRUE , signif.light, insignif.light), width=0) +
            geom_errorbar(aes(ymin=co - se, ymax=co + se), lwd=2.5, 
                          colour = ifelse(signif == TRUE, signif.medium, insignif.medium), width=0) +
            geom_point(size = 3, pch = ifelse(signif == TRUE, 21, 22), fill = ifelse(signif == TRUE, signif.dark, insignif.dark)) +
            coord_flip() + 
            theme_bw() 
        if (length(models) == 1) {
            #cat(length(models))
            p <- p +   ggtitle("Model")
        } else if (length(models) > 1) {
            p <- p +  facet_wrap(~co.names, strip.position="left", nrow=length(dataframe),scales = "free_y") +
                ggtitle("Models")
        }
        
    }
    
    print(p)
    #meassages to paste
    # add messages to p and print at the end
    if (use.se[i] == TRUE && length(models[[i]]@se) > 0) {
        note <- "Bars denote SEs."
        message(paste0("Model ", i, 
                       ": bars denote one (inner) resp. two (outer) standard errors."))
    } else if (length(models[[i]]@ci.low) == 0 && length(models[[i]]@se) > 0) {
        note <- "Bars denote CIs."
        message(paste0("Model ", i, ": bars denote 0.5 (inner) resp. ", ci.level, 
                       " (outer) confidence intervals (computed from standard errors)."))
    } else {
        note <- "Bars denote CIs."
        message(paste0("Model ", i, ": bars denote  ", ci.level, 
                       " confidence intervals."))
    }
    
    if (!is.null(custom.note)) {
        note <- custom.note
    }
    
    if (length(co) == 0) {
        stop(paste("No coefficients available. Was the 'omit.coef' argument", 
                   "misspecified? If coefficients were renamed using the",
                   "'custom.coef.names' argument, 'omit.coef' refers to the renamed",
                   "coefficients."))
    }
    
    if (!is.null(file) && !is.na(file)) {
        dev.off()
    }
}

