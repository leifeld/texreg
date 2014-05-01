
# function which plots coefficients and confidence intervals ('coefplot')
coefplot <- function(labels, estimates, lower.inner = NULL, 
    upper.inner = NULL, lower.outer = NULL, upper.outer = NULL, 
    signif.outer = TRUE, xlab = "Coefficients and confidence intervals", 
    main = "Coefficient plot", vertical.lines = TRUE, cex = 2.5, 
    lwd.inner = 7, lwd.outer = 5, signif.light = "#fbc9b9", 
    signif.medium = "#f7523a", signif.dark = "#bd0017", 
    insignif.light = "#c5dbe9", insignif.medium = "#5a9ecc", 
    insignif.dark = "#1c5ba6", ...) {
  
  # check consistency of arguments
  if (length(table(c(length(estimates), length(lower.outer), 
      length(upper.outer), length(labels)))) > 1) {
    stop("Vectors should have the same length.")
  }
  
  # define shortcut variables outer, mn, mx, steps, and num
  if (!is.null(lower.outer) && !is.null(upper.outer)) {
    outer <- TRUE
  } else {
    outer <- FALSE
  }
  if (!is.null(lower.inner) && !is.null(upper.inner)) {
    inner <- TRUE
  } else {
    inner <- FALSE
  }
  if (outer == TRUE) {
    mn <- min(lower.outer)
    mx <- max(upper.outer)
  } else if (inner == TRUE) {
    mn <- min(lower.inner)
    mx <- max(upper.inner)
  } else {
    stop("Either inner or outer bounds must be provided.")
  }
  if (mn > 0) {
    mn <- 0
  }
  if (mx < 0) {
    mx <- 0
  }
  steps <- floor(mn):ceiling(mx)
  num <- length(labels)
  
  # which terms are significant?
  if (class(signif.outer) == "logical" && 
      length(signif.outer) == length(estimates)) {
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
  
  # create plot; compute left margin; add labels to left margin
  inches.to.lines <- (par("mar") / par("mai"))[1]
  lab.width <- max(strwidth(labels, units = "inches") * inches.to.lines)
  
  par(mar = c(5, 1 + lab.width, 4, 2) + 0.1)
  plot(0, 0, xlim = c(mn, mx), ylim = c(1, length(labels)), xlab = xlab, 
      main = main, ylab = "", axes = FALSE, type = "n", ...)
  axis(side = 2, at = 1:num, labels = FALSE, las = 2, 
      tck = 0, lty = 0)
  if (any(signif)) {
    mtext(side = 2, text = labels[which(signif)], at = num + 1 - which(signif), 
        col = signif.dark, las = 2)
  }
  if (any(!signif)) {
    mtext(side = 2, text = labels[which(!signif)], at = num + 1 - which(!signif), 
        col = insignif.dark, las = 2)
  }
  axis(side = 1, las = 0, lty = 0)
  
  # add vertical gray lines
  if (vertical.lines == TRUE) {
    zeros <- rep(0, length(steps))
    ends <- rep((num + 1), length(steps))
    segments(steps, zeros, steps, ends, col = "gray", lwd = 1)
  }
  segments(0, 0, 0, num + 1, lwd = 4, col = "grey75")
  
  # draw outer CIs
  if (outer == TRUE) {
    if (any(signif)) {
      segments(lower.outer[signif], num + 1 - which(signif), 
          upper.outer[signif], num + 1 - which(signif), lwd = lwd.outer, 
          col = signif.light)
    }
    if (any(!signif)) {
      segments(lower.outer[!signif], num + 1 - which(!signif), 
          upper.outer[!signif], num + 1 - which(!signif), lwd = lwd.outer, 
          col = insignif.light)
    }
  }
  
  # draw inner CIs
  if (inner == TRUE) {
    if (any(signif)) {
      segments(lower.inner[signif], num + 1 - which(signif), 
          upper.inner[signif], num + 1 - which(signif), lwd = lwd.inner, 
          col = signif.medium)
    }
    if (any(!signif)) {
      segments(lower.inner[!signif], num + 1 - which(!signif), 
          upper.inner[!signif], num + 1 - which(!signif), lwd = lwd.inner, 
          col = insignif.medium)
    }
    if (outer == FALSE) {

    }
  }
  
  # draw point estimates
  if (any(signif)) {
    points(estimates[signif], num + 1 - which(signif), col = signif.dark, 
        pch = 20, cex = cex)
  }
  if (any(!signif)) {
    points(estimates[!signif], num + 1 - which(!signif), col = insignif.dark, 
        pch = 20, cex = cex)
  }
}


# plotreg function
plotreg <- function(l, file = NA, custom.model.names = NULL, 
    custom.coef.names = NULL, custom.note = NULL, override.coef = 0, 
    override.se = 0, override.pval = 0, omit.coef = NA, reorder.coef = NULL, 
    ci.level = 0.95, use.se = FALSE, mfrow = TRUE, vertical.lines = TRUE, 
    cex = 2.5, lwd.inner = 7, lwd.outer = 5, signif.light = "#fbc9b9", 
    signif.medium = "#f7523a", signif.dark = "#bd0017", 
    insignif.light = "#c5dbe9", insignif.medium = "#5a9ecc", 
    insignif.dark = "#1c5ba6", ...) {
  
  if (!is.na(omit.coef) && !is.character(omit.coef)) {
    stop("omit.coef must be a character string!")
  }
  
  if (length(use.se) == 1) {
    use.se <- rep(use.se, length(l))
  }
  
  # extract texreg objects and override data
  models <- get.data(l, ...)
  models <- override(models, override.coef, override.se, override.pval)
  
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
  if (!is.na(file)) {
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
    remove.rows <- grep(omit.coef, models[[i]]@coef.names)
    if (is.na(omit.coef)) {
      remove.rows <- length(models[[i]]@coef) + 1
    }
    
    # aggregate confidence intervals or SEs
    if (use.se[i] == TRUE && length(models[[i]]@se) > 0) {
      co <- models[[i]]@coef
      co.names <- models[[i]]@coef.names
      se <- models[[i]]@se
      pv <- models[[i]]@pvalues
      co <- co[-remove.rows]
      co.names <- co.names[-remove.rows]
      se <- se[-remove.rows]
      lower.inner <- co - se
      upper.inner <- co + se
      lower.outer <- co - 2 * se
      upper.outer <- co + 2 * se
      if (length(pv) > 0) {
        pv <- pv[-remove.rows]
        signif.outer <- pv < (1 - ci.level)
        
        # sort according to reorder.coef argument
        dataframe <- data.frame(co.names, co, se, pv)
        dataframe <- reorder(dataframe, reorder.coef)
        pv <- dataframe[, 4]
      } else {
        # sort according to reorder.coef argument
        dataframe <- data.frame(co.names, co, se)
        dataframe <- reorder(dataframe, reorder.coef)
      }
      co.names <- dataframe[, 1]
      co <- dataframe[, 2]
      se <- dataframe[, 3]
      note <- "Bars denote SEs."
      cat(paste0("Model ", i, 
        ": bars denote one (inner) resp. two (outer) standard errors.\n"))
    } else if (length(models[[i]]@ci.low) == 0 && length(models[[i]]@se) > 0) {
      co <- models[[i]]@coef
      co.names <- models[[i]]@coef.names
      se <- models[[i]]@se
      co <- co[-remove.rows]
      co.names <- co.names[-remove.rows]
      se <- se[-remove.rows]
      z.inner <- qnorm(1 - ((1 - 0.5) / 2))
      lower.inner <- co - (z.inner * se)
      upper.inner <- co + (z.inner * se)
      z.outer <- qnorm(1 - ((1 - ci.level) / 2))
      lower.outer <- co - (z.outer * se)
      upper.outer <- co + (z.outer * se)
      
      # sort according to reorder.coef argument
      dataframe <- data.frame(co.names, co, lower.inner, upper.inner, 
          lower.outer, upper.outer)
      dataframe <- reorder(dataframe, reorder.coef)
      co.names <- dataframe[, 1]
      co <- dataframe[, 2]
      lower.inner <- dataframe[, 3]
      upper.inner <- dataframe[, 4]
      lower.outer <- dataframe[, 5]
      upper.outer <- dataframe[, 6]
      
      signif.outer <- TRUE
      note <- "Bars denote CIs."
      cat(paste0("Model ", i, ": bars denote 0.5 (inner) resp. ", ci.level, 
          " (outer) confidence intervals (computed from standard errors).\n"))
    } else {
      co <- models[[i]]@coef
      co.names <- models[[i]]@coef.names
      co <- co[-remove.rows]
      co.names <- co.names[-remove.rows]
      lower.outer <- models[[i]]@ci.low
      upper.outer <- models[[i]]@ci.up
      lower.inner <- NULL
      upper.inner <- NULL
      lower.outer <- lower.outer[-remove.rows]
      upper.outer <- upper.outer[-remove.rows]
      
      # sort according to reorder.coef argument
      dataframe <- data.frame(co.names, co, lower.outer, upper.outer)
      dataframe <- reorder(dataframe, reorder.coef)
      co.names <- dataframe[, 1]
      co <- dataframe[, 2]
      lower.outer <- dataframe[, 3]
      upper.outer <- dataframe[, 4]
      
      signif.outer <- TRUE
      note <- "Bars denote CIs."
      cat(paste0("Model ", i, ": bars denote 0.5 (inner) resp. ", ci.level, 
          " (outer) confidence intervals.\n"))
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
    
    # plot
    coefplot(
        labels = co.names, 
        estimates = co,
        lower.inner = lower.inner,
        upper.inner = upper.inner,
        lower.outer = lower.outer,
        upper.outer = upper.outer,
        signif.outer = signif.outer,
        xlab = note, 
        main = model.names[i], 
        vertical.lines = vertical.lines, 
        cex = cex, 
        lwd.inner = lwd.inner, 
        lwd.outer = lwd.outer, 
        signif.light = signif.light,
        signif.medium = signif.medium, 
        signif.dark = signif.dark, 
        insignif.light = insignif.light, 
        insignif.medium = insignif.medium, 
        insignif.dark = insignif.dark, 
        ...
    )
  }
  
  if (!is.na(file)) {
    dev.off()
  }
}

