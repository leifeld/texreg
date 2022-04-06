context("JSS article 2013")
suppressPackageStartupMessages(library("texreg"))

# example models from ?lm
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
weight <- c(ctl, trt)
m1 <- lm(weight ~ group)
m2 <- lm(weight ~ group - 1)

test_that("texreg returns output as in the JSS 2013 article", {

  # Simple screenreg example
  expect_equal(output <- screenreg(list(m1, m2)),
               readRDS("../files/jss_screenreg_lm.RDS"))
  # saveRDS(output, "../files/jss_screenreg_lm.RDS")

  # texreg example with dcolumn and booktabs usage and table float options
  expect_equal(output <- texreg(list(m1, m2),
                      dcolumn = TRUE,
                      booktabs = TRUE,
                      use.packages = FALSE,
                      label = "tab:3",
                      caption = "Two linear models.",
                      float.pos = "bh"),
               readRDS("../files/jss_texreg_dcolumn_booktabs.RDS"))
  # saveRDS(output, "../files/jss_texreg_dcolumn_booktabs.RDS")

  # Bold coefficients, custom note, omit.coef, and coefficient customization
  # (difference to JSS: dollar signs around GOF values; but appearance otherwise identical)
  expect_equal(output <- texreg(list(m1, m2),
                                label = "tab:4",
                                caption = "Bolded coefficients, custom notes, three digits.",
                                float.pos = "h",
                                return.string = TRUE,
                                bold = 0.05,
                                stars = 0,
                                custom.note = "Coefficients with $p < 0.05$ in \\textbf{bold}.",
                                digits = 3,
                                leading.zero = FALSE,
                                omit.coef = "Inter"),
               readRDS("../files/jss_texreg_bold_customnote_digits.RDS"))
  # saveRDS(output, "../files/jss_texreg_bold_customnote_digits.RDS")

  # GLS example; custom names, reordering, single.row, 'extract' arguments
  # (difference to JSS: the paper reports results using 'no.margin = TRUE', but it's not in the code example)
  # (difference to JSS: the version used in the paper counts 11 places left of the right bracket; this is now correctly counted as 9)
  expect_equal({
      library("nlme")
      m3 <- gls(follicles ~ sin(2 * pi * Time) + cos(2 * pi * Time), Ovary,
                correlation = corAR1(form = ~ 1 | Mare))

      table <- texreg(list(m1, m3),
                      custom.coef.names = c(
                        "Intercept",
                        "Control",
                        "$\\sin(2 \\cdot \\pi \\cdot \\mbox{time})$",
                        "$\\cos(2 \\cdot \\pi \\cdot \\mbox{time})$"
                      ),
                      custom.model.names = c("OLS model", "GLS model"),
                      reorder.coef = c(1, 3, 4, 2),
                      caption = "Multiple model types, custom names, and single row.",
                      label = "tab:5",
                      stars = c(0.01, 0.001),
                      dcolumn = TRUE,
                      booktabs = TRUE,
                      use.packages = FALSE,
                      no.margin = TRUE,
                      single.row = TRUE,
                      include.adjrs = FALSE,
                      include.bic = FALSE)
    },
    readRDS("../files/jss_texreg_gls.RDS")
  )
  # saveRDS(table, "../files/jss_texreg_gls.RDS")

  # How to use "robust" standard errors with texreg
  expect_equal({
      library("sandwich")
      library("lmtest")
      hc <- vcovHC(m2)
      ct <- coeftest(m2, vcov = hc)
      se <- ct[, 2]
      pval <- ct[, 4]
      output <- texreg(m2, override.se = se, override.pvalues = pval)
    },
    readRDS("../files/jss_texreg_robust.RDS")
  )
  # saveRDS(output, "../files/jss_texreg_robust.RDS")

  # Creating Word-readable HTML files using htmlreg
  expect_equal({
      output <- htmlreg(list(m1, m2, m3),
                        inline.css = FALSE,
                        doctype = TRUE,
                        html.tag = TRUE,
                        head.tag = TRUE,
                        body.tag = TRUE)
    },
    readRDS("../files/jss_htmlreg_word.RDS")
  )
  # saveRDS(output, "../files/jss_htmlreg_word.RDS")

  # Compatibility with Markdown
  expect_equal({
      output <- htmlreg(list(m1, m2, m3), star.symbol = "\\*", center = TRUE)
    },
    readRDS("../files/jss_htmlreg_markdown.RDS")
  )
  # saveRDS(output, "../files/jss_htmlreg_markdown.RDS")

  # How to write a complete extension for linear models
  expect_equal({
      extract.lm <- function(model, include.rsquared = TRUE,
                             include.adjrs = TRUE, include.nobs = TRUE, ...) {

        s <- summary(model, ...)
        names <- rownames(s$coef)
        co <- s$coef[, 1]
        se <- s$coef[, 2]
        pval <- s$coef[, 4]

        gof <- numeric()
        gof.names <- character()
        gof.decimal <- logical()
        if (include.rsquared == TRUE) {
          rs <- s$r.squared
          gof <- c(gof, rs)
          gof.names <- c(gof.names, "R$^2$")
          gof.decimal <- c(gof.decimal, TRUE)
        }
        if (include.adjrs == TRUE) {
          adj <- s$adj.r.squared
          gof <- c(gof, adj)
          gof.names <- c(gof.names, "Adj.\\ R$^2$")
          gof.decimal <- c(gof.decimal, TRUE)
        }
        if (include.nobs == TRUE) {
          n <- nobs(model)
          gof <- c(gof, n)
          gof.names <- c(gof.names, "Num.\\ obs.")
          gof.decimal <- c(gof.decimal, FALSE)
        }

        tr <- createTexreg(
          coef.names = names,
          coef = co,
          se = se,
          pvalues = pval,
          gof.names = gof.names,
          gof = gof,
          gof.decimal = gof.decimal
        )
        return(tr)
      }
      setMethod("extract", signature = className("lm", "stats"), definition = extract.lm)
    },
    "extract"
  )
})
