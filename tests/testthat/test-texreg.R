context("texreg family functions")
library("texreg")

data("iris")
model1 <- lm(Sepal.Width ~ Petal.Width, data = iris)

test_that("no stars, single.row, and dcolumn work together", {
  expect_match(texreg(model1,
                      single.row = TRUE,
                      dcolumn = TRUE,
                      stars = numeric()),
               regexp = "l D{\\)}{\\)}{11\\)0}}",
               perl = TRUE)
})

test_that("sanity checks and arguments work in texreg function", {
  expect_warning(texreg(model1, dcolumn = TRUE, bold = 0.05),
                 "cannot be used at the same time. Switching off 'dcolumn'. You should also consider setting stars = 0.")
  expect_warning(texreg(model1, dcolumn = TRUE, bold = 0.05, stars = 0), "cannot be used at the same time. Switching off 'dcolumn'.$")
  expect_warning(texreg(model1, dcolumn = TRUE, bold = 0.05, stars = numeric()), "cannot be used at the same time. Switching off 'dcolumn'.$")
  expect_warning(texreg(model1, longtable = TRUE, sideways = TRUE), "Switching off 'sideways'.$")
  expect_match(texreg(model1, longtable = TRUE), "\\\\usepackage{longtable}\n\n", perl = TRUE)
  expect_match(texreg(model1, longtable = TRUE, float.pos = "tbhp"), "\\\\begin{longtable}\\[tbhp\\]{l c}\n", perl = TRUE)
  expect_match(texreg(model1, scalebox = 0.5), "\\\\usepackage{graphicx}\n\n", perl = TRUE)
  expect_match(texreg(model1, scalebox = 0.5), "\\\\scalebox{0.5}{", perl = TRUE)
  expect_warning(texreg(model1, longtable = TRUE, scalebox = 0.5),
                 "'longtable' and 'scalebox' are not compatible. Setting scalebox = NULL.$")
  expect_match(texreg(model1, sideways = TRUE), "\\\\usepackage{rotating}\n\n", perl = TRUE)
  expect_match(texreg(model1, sideways = TRUE), "\\\\begin{sidewaystable}\n", perl = TRUE)
  expect_match(texreg(model1, booktabs = TRUE), "\\\\usepackage{booktabs}\n\n", perl = TRUE)
  expect_match(texreg(model1, booktabs = TRUE), "\\\\toprule\n", perl = TRUE)
  expect_match(texreg(model1, booktabs = TRUE), "\\\\midrule\n", perl = TRUE)
  expect_match(texreg(model1, booktabs = TRUE), "\\\\bottomrule\n", perl = TRUE)
  expect_match(texreg(model1, ci.force = TRUE), "& \\$ \\[ 3.19;  3.43\\]\\$ \\\\\\\\\n", perl = TRUE)
  expect_match(texreg(model1, ci.force = TRUE, single.row = TRUE),
               "\\(Intercept\\) & \\$3.31 \\\\; \\[ 3\\.19;  3\\.43\\]\\^{.}\\$  \\\\\\\\\n",
               perl = TRUE)
  expect_match(texreg(model1, center = FALSE), "\\\\end{table}$", perl = TRUE)
  expect_match(texreg(model1, float.pos = ""), "\\\\begin{table}\n", perl = TRUE)
  expect_match(texreg(model1, float.pos = "tb"), "\\\\begin{table}\\[tb\\]\n", perl = TRUE)
  expect_match(texreg(model1, fontsize = "footnotesize"), "\\\\begin{footnotesize}\n", perl = TRUE)
  expect_match(texreg(model1, fontsize = "footnotesize"), "\\\\end{footnotesize}\n", perl = TRUE)
  expect_match(texreg(model1, fontsize = "footnotesize"), "\\\\multicolumn\\{2\\}{l}{\\\\tiny{", perl = TRUE)
  expect_match(texreg(model1, fontsize = "normalsize"), "\\\\multicolumn\\{2\\}{l}{\\\\scriptsize{", perl = TRUE)
  expect_match(texreg(model1, fontsize = "large"), "\\\\multicolumn\\{2\\}{l}{\\\\footnotesize{", perl = TRUE)
  expect_match(texreg(model1, fontsize = "Large"), "\\\\multicolumn\\{2\\}{l}{\\\\small{", perl = TRUE)
  expect_match(texreg(model1, fontsize = "LARGE"), "\\\\multicolumn\\{2\\}{l}{\\\\normalsize{", perl = TRUE)
  expect_match(texreg(model1, fontsize = "huge"), "\\\\multicolumn\\{2\\}{l}{\\\\large{", perl = TRUE)
  expect_match(texreg(model1, fontsize = "Huge"), "\\\\multicolumn\\{2\\}{l}{\\\\Large{", perl = TRUE)
  expect_match(texreg(model1, fontsize = "footnotesize", longtable = TRUE), "\\\\begin{footnotesize}\n\\\\begin{longtable}", perl = TRUE)
  expect_match(texreg(model1, caption.above = TRUE), "\\\\begin{table}\n\\\\caption{Statistical models}\n", perl = TRUE)
  expect_match(texreg(model1, custom.note = "my new note"), "\\\\multicolumn\\{2\\}{l}{\\\\scriptsize{my new note}}", perl = TRUE)
  expect_match(texreg(model1, custom.note = ""), "\\\\hline\n\\\\end{tabular}\n", perl = TRUE)
  expect_match(texreg(model1, custom.note = "my new note", longtable = TRUE),
               "\\\\multicolumn\\{2\\}{l}{\\\\scriptsize{my new note}}\\\\\\\\\n",
               perl = TRUE)
  expect_match(texreg(model1, longtable = TRUE, caption.above = TRUE),
               "\\\\begin{longtable}{l c}\n\\\\caption{Statistical models}\n",
               perl = TRUE)
  expect_gt({
      suppressMessages(texreg(model1, file = "../files/temp.txt"))
      fs <- file.size("../files/temp.txt")
      unlink("../files/temp.txt")
      fs
    }, 400)
  expect_error(texreg(model1, file = 12), "The 'file' argument must be a character string.")
  expect_match(texreg(model1,
                      custom.columns = list("column 1" = c("yes", "no"),
                                            c(11, 13)),
                      custom.col.pos = c(1, 2)),
               "column 1 &  &  & Model 1 \\\\\\\\\n",
               perl = TRUE)
  expect_match(texreg(list(model1, model1),
                      custom.columns = list("column 1" = c("yes", "no"),
                                            c(11, 13)),
                      custom.col.pos = c(1, 3)),
               "column 1 &  & Model 1 &  & Model 2 \\\\\\\\\n",
               perl = TRUE)
  expect_match(texreg(list(model1, model1),
                      custom.columns = list("column 1" = c("yes", "no"),
                                            c(11, 13)),
                      custom.col.pos = c(1, 3)),
               "yes & \\(Intercept\\) & \\$3.31\\^{.{3}}\\$  & 11 & \\$3\\.31\\^{.{3}}\\$  \\\\\\\\\n",
               perl = TRUE)
  expect_match(texreg(model1, dcolumn = TRUE, custom.columns = list("column 1" = 1:2)),
               " & column 1 & \\\\multicolumn\\{1\\}{c}{Model 1} \\\\\\\\\n",
               perl = TRUE)
  expect_match(capture.output(print(texreg(model1)))[22], "\\\\end\\{table\\}")

  # custom.gof.rows with correct dcolumn width and without dollar signs for text
  expect_match(texreg(model1, dcolumn = TRUE, custom.gof.rows = list("Fixed effects" = 12345)), "Fixed effects \\& 12345       \\\\")
  expect_match(texreg(model1, dcolumn = TRUE, custom.gof.rows = list("Fixed effects" = 12345)), "\\\\begin[{]tabular[}][{]l D[{]\\.[}][{]\\.[}][{]5\\.5[}][}]")
  expect_match(texreg(model1, dcolumn = FALSE, custom.gof.rows = list("Fixed effects" = 12345)), "Fixed effects \\& \\$12345\\$       \\\\")
  expect_match(texreg(model1, dcolumn = FALSE, custom.gof.rows = list("Fixed effects" = TRUE)), "Fixed effects \\& TRUE          \\\\")
  expect_match(texreg(model1, dcolumn = FALSE, custom.gof.rows = list("Fixed effects" = "yes")), "Fixed effects \\& yes           \\\\")
  expect_match(texreg(model1, dcolumn = TRUE, custom.gof.rows = list("Fixed effects" = TRUE)), "Fixed effects \\& \\\\multicolumn\\{1\\}\\{c\\}\\{TRUE\\} \\\\")
  expect_match(texreg(model1, dcolumn = TRUE, custom.gof.rows = list("Fixed effects" = TRUE)), "3\\.5")
  expect_match(texreg(model1, dcolumn = TRUE, custom.gof.rows = list("Fixed effects" = "yes")), "Fixed effects \\& \\\\multicolumn\\{1\\}\\{c\\}\\{yes\\} \\\\")
  expect_match(texreg(model1, dcolumn = TRUE, custom.gof.rows = list("Fixed effects" = "yes")), "3\\.5")
  expect_match(texreg(model1, dcolumn = TRUE, custom.gof.rows = list("Fixed effects" = 1234567)), "7\\.5")

  # groups
  expect_match(texreg(model1, groups = list("First group" = 1, "Second group" = 2), single.row = FALSE),
               "Second group      \\&               \\\\")
  expect_match(texreg(model1, groups = list("First group" = 1, "Second group" = 2), single.row = FALSE),
               "\n\\\\quad Petal\\.Width")
})

test_that("arguments work in screenreg function", {
  expect_error(screenreg(model1, outer.rule = TRUE), "outer.rule must be a character.")
  expect_error(screenreg(model1, outer.rule = "=="), "outer.rule must be a character of maximum length 1.")
  expect_match(screenreg(model1, outer.rule = ""), "^\n     ")
  expect_error(screenreg(model1, inner.rule = TRUE), "inner.rule must be a character.")
  expect_error(screenreg(model1, inner.rule = "--"), "inner.rule must be a character of maximum length 1.")
  expect_match(screenreg(model1, inner.rule = ""), "Model 1   \n\\(Intercept")
  expect_match(screenreg(model1, custom.note = ""), "==\n$")
  expect_match(screenreg(model1, custom.note = "my note"), "==\nmy note\n$")
  expect_match(screenreg(model1, custom.note = "%stars. my note"), "0.05. my note\n$")
  expect_error(screenreg(model1, file = TRUE), "The 'file' argument must be a character string.$")
  expect_gt({
    suppressMessages(screenreg(model1, file = "../files/temp.txt"))
    fs <- file.size("../files/temp.txt")
    unlink("../files/temp.txt")
    fs
  }, 300)
  expect_match(screenreg(list(model1, model1),
                         custom.gof.rows = list(test = c("yes", "no"),
                                                "test 2" = 4:5,
                                                "test 3" = c(2.4, 3.7))),
                         "--\ntest\\s+yes\\s+no\\s+\ntest 2\\s+4\\s+5\\s+\ntest 3\\s+2.40\\s+3.70\\s+\nR\\^2")

  # groups
  expect_match(screenreg(model1, groups = list("First group" = 1, "Second group" = 2), single.row = FALSE),
               "\\nFirst group                \\n")
  expect_match(screenreg(model1, groups = list("First group" = 1, "Second group" = 2), single.row = FALSE),
               "\\n    Petal\\.Width")
})

test_that("knitreg function works", {
  with_mock(requireNamespace = function (package, ...) {
    ifelse(package == "knitr", return(FALSE), return(TRUE))
  }, {
    expect_error(knitreg(list(model1, model1)), regex = "knitreg requires the 'knitr' package to be installed")
  })
  with_mock(requireNamespace = function (package, ...) {
      ifelse(package == "rmarkdown", return(FALSE), return(TRUE))
    }, {
    expect_error(knitreg(list(model1, model1)), regex = "knitreg requires the 'rmarkdown' package to be installed")
  })
  skip_if_not_installed("knitr", minimum_version = "1.22")
  require("knitr")
  skip_if_not_installed("rmarkdown", minimum_version = "1.12")
  require("rmarkdown")
  expect_match(knitreg(model1), "Petal.Width")
})

test_that("matrixreg function works", {
  expect_equal(nrow(matrixreg(model1)), 8)
  expect_equal(ncol(matrixreg(model1)), 2)
  expect_equal(class(matrixreg(model1))[1], "matrix")
  expect_error(matrixreg(model1, custom.coef.map = list("Intercept" = "My intercept")), "None of the coefficient")
  mr <- matrixreg(model1, custom.coef.map = list("(Intercept)" = "My intercept"))
  expect_equal(nrow(mr), 6)
  expect_match(mr[2, 1], "My intercept")
})

test_that("htmlreg function works", {
  # groups
  expect_match(htmlreg(model1, groups = list("First group" = 1, "Second group" = 2), single.row = FALSE),
               "\\&nbsp;\\&nbsp;\\&nbsp;\\&nbsp;\\&nbsp;Petal\\.Width")
  expect_match(htmlreg(model1, groups = list("First group" = 1, "Second group" = 2), single.row = FALSE),
               "Second group")
})