context("Test of texreg et al. function arguments")
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
})