context("texreg family functions")
suppressPackageStartupMessages(library("texreg"))

data("iris")
model1 <- lm(Sepal.Width ~ Petal.Width, data = iris)

test_that("custom.header argument works", {
  # screenreg
  sr1 <- screenreg(list(model1, model1, model1, model1, model1, model1),
                   custom.header = list("ab" = 1:2, "def" = 4:6))
  expect_match(sr1, "\n                       ab                                        def               \n             ----------------------              ----------------------------------\n             Model 1")
  sr2 <- screenreg(list(model1, model1, model1, model1, model1, model1),
                   custom.header = list("ab" = 1:2, "def" = 4:6),
                   custom.columns = list("first" = c("one", "two"), "second" = c("three", "four")),
                   custom.col.pos = c(1, 6))
  expect_match(sr2, " ---------------------- ")
  expect_match(sr2, " ------------------------------------------\n")
  expect_match(sr2, "-              -")
  expect_match(sr2, "\nfirst               Model 1     Model 2     Model 3     Model 4     second  Model 5     Model 6   \n")
  sr3 <- screenreg(list(model1, model1, model1),
                   custom.header = list("ab" = 1:2, "defghijklmnopqrstuvxwyz" = 3),
                   single.row = TRUE, groups = list("Group 1" = 1, "Group 2" = 2))
  expect_match(sr3, " defghijklmnopqrst\n")
  expect_match(sr3, "\n    Petal.Width")

  # texreg
  tr1 <- texreg(list(model1, model1, model1, model1, model1, model1),
                custom.header = list("ab" = 1:2, "def" = 4:6))
  expect_match(tr1, "\\n \\& \\\\multicolumn\\{2\\}\\{c\\}\\{ab\\} \\& \\& \\\\multicolumn\\{3\\}\\{c\\}\\{def\\} \\\\\\\\\\n\\\\cline\\{2-3\\} \\\\cline\\{5-7\\}\\n")
  tr2 <- texreg(list(model1, model1, model1, model1, model1, model1),
                custom.header = list("ab" = 1:2, "def" = 4:6),
                custom.columns = list("first" = c("one", "two"), "second" = c("three", "four")),
                custom.col.pos = c(1, 6))
  expect_match(tr2, "\\n & & \\\\multicolumn\\{2\\}\\{c\\}\\{ab\\} \\& \\& \\\\multicolumn\\{4\\}\\{c\\}\\{def\\} \\\\\\\\\\n\\\\cline\\{3-4\\} \\\\cline\\{6-9\\}\\n")
  tr3 <- texreg(list(model1, model1, model1),
                custom.header = list("ab" = 1:2, "defghijklmnopqrstuvxwyz" = 3),
                single.row = TRUE, groups = list("Group 1" = 1, "Group 2" = 2),
                booktabs = TRUE, longtable = TRUE, siunitx = TRUE)
  expect_match(tr3, "\\n\\\\toprule\\n \\& \\\\multicolumn\\{2\\}\\{c\\}\\{ab\\} \\& \\\\multicolumn\\{1\\}\\{c\\}\\{defghijklmnopqrstuvxwyz\\} \\\\\\\\\\n\\\\cmidrule\\(lr\\)\\{2-3\\} \\\\cmidrule\\(lr\\)\\{4-4\\}\\n")
  expect_error(texreg(list(model1, model1), custom.header = c("a", "b")), "must be a named list of numeric vectors")
  expect_error(texreg(list(model1, model1), custom.header = list("ab" = c(1, NA))), "NA values are not permitted in 'custom.header'")
  expect_error(texreg(list(model1, model1), custom.header = list("ab" = c(1, 1.5))), "must be provided as integer values")
  expect_error(texreg(list(model1, model1), custom.header = list("ab" = c(1, 3))), "must be between 1 and the number of models")
  expect_error(texreg(list(model1, model1), custom.header = list("ab" = c(2, 1))), "must be strictly increasing")
  expect_error(texreg(list(model1, model1, model1), custom.header = list("ab" = c(1, 3))), "must have strictly consecutive column indices")

  # htmlreg
  hr1 <- htmlreg(list(model1, model1, model1, model1, model1, model1),
                 custom.header = list("ab" = 1:2, "def" = 4:6),
                 inline.css = FALSE)
  expect_match(hr1, "<tr>\\n<th>\\&nbsp;<\\/th>\\n<th class = \"rule\" colspan=\"2\">ab<\\/th>\\n<th>\\&nbsp;<\\/th>\\n<th class = \"rule\" colspan=\"3\">def<\\/th>\\n<\\/tr>")
  hr2 <- htmlreg(list(model1, model1, model1, model1, model1, model1),
                 custom.header = list("ab" = 1:2, "def" = 4:6),
                 custom.columns = list("first" = c("one", "two"), "second" = c("three", "four")),
                 custom.col.pos = c(1, 6), inline.css = FALSE, body.tag = TRUE)
  expect_match(hr2, "<tr>\\n<th>\\&nbsp;<\\/th>\\n<th>\\&nbsp;<\\/th>\\n<th class = \"rule\" colspan=\"2\">ab<\\/th>\\n<th>\\&nbsp;<\\/th>\\n<th class = \"rule\" colspan=\"4\">def<\\/th>\\n<\\/tr>")
})

test_that("ci.test argument works", {
  mr <- matrixreg(list(model1, model1), ci.force = TRUE, ci.test = 3.3)
  expect_match(mr[2, 2], "  3.31        ")
  expect_match(mr[2, 3], "  3.31        ")
  expect_match(mr[4, 2], " -0.21 \\*      ")
  expect_match(mr[4, 3], " -0.21 \\*      ")
  mr <- matrixreg(list(model1, model1), ci.force = TRUE, ci.test = 0)
  expect_match(mr[2, 2], "  3.31 \\*      ")
  expect_match(mr[2, 3], "  3.31 \\*      ")
  expect_match(mr[4, 2], " -0.21 \\*      ")
  expect_match(mr[4, 3], " -0.21 \\*      ")
  mr <- matrixreg(list(model1, model1), ci.force = TRUE, ci.test = NA)
  expect_match(mr[2, 2], "  3.31        ")
  expect_match(mr[2, 3], "  3.31        ")
  expect_match(mr[4, 2], " -0.21        ")
  expect_match(mr[4, 3], " -0.21        ")
  mr <- matrixreg(list(model1, model1), ci.force = TRUE, ci.test = c(0, 3.3))
  expect_match(mr[2, 2], "  3.31 \\*      ")
  expect_match(mr[2, 3], "  3.31        ")
  expect_match(mr[4, 2], " -0.21 \\*      ")
  expect_match(mr[4, 3], " -0.21 \\*      ")
  mr <- matrixreg(list(model1, model1), ci.force = TRUE, ci.test = c(NA, 3.3))
  expect_match(mr[2, 2], "  3.31        ")
  expect_match(mr[2, 3], "  3.31        ")
  expect_match(mr[4, 2], " -0.21        ")
  expect_match(mr[4, 3], " -0.21 \\*      ")
  mr <- matrixreg(list(model1, model1), ci.force = TRUE, ci.test = c(NA, NA))
  expect_match(mr[2, 2], "  3.31        ")
  expect_match(mr[2, 3], "  3.31        ")
  expect_match(mr[4, 2], " -0.21        ")
  expect_match(mr[4, 3], " -0.21        ")

  # check multiple NA values
  sr <- screenreg(list(model1, model1),
                  custom.model.names = c("Model 1", "Model 2"),
                  ci.force = TRUE,
                  ci.test = c(NA, NA)
  )
  expect_false(grepl("[*]", sr))
  sr2 <- screenreg(list(model1, model1),
                   custom.model.names = c("Model 1", "Model 2"),
                   ci.force = TRUE,
                   ci.test = c(NA, 0)
  )
  expect_equal(attr(regexpr("[*]", sr2), "match.length"), 1)
})

test_that("threeparttable and custom.note arguments work in the texreg function", {
  tr <- texreg(model1, threeparttable = TRUE)
  expect_match(tr, "\\usepackage\\{threeparttable\\}\\n\\n")
  expect_match(tr, "\\\\begin\\{center\\}\\n\\\\begin\\{threeparttable\\}\\n\\\\begin\\{tabular\\}", perl = TRUE)
  expect_match(tr, "\\\\end{tabular}\\n\\\\begin\\{tablenotes\\}\\[flushleft\\]\\n\\\\scriptsize\\{.+\\}\\n\\\\end\\{tablenotes\\}\\n\\\\end\\{threeparttable\\}\\n\\\\caption", perl = TRUE)
  tr <- texreg(model1, threeparttable = TRUE, custom.note = "\n\\item %stars.\\\\\n\\item Second item.\n")
  expect_match(tr, "\\\\end{tabular}\\n\\\\begin\\{tablenotes\\}\\[flushleft\\]\\n\\\\scriptsize\\{\\n.+\\n", perl = TRUE)
  expect_match(tr, "\\\\\\n\\\\item Second item\\.\\n\\}\\n\\\\end\\{tablenotes\\}", perl = TRUE)
  tr <- texreg(model1, threeparttable = TRUE, sideways = TRUE)
  expect_match(tr, "\\\\begin\\{sidewaystable\\}\\n\\\\begin\\{center\\}\\n\\\\begin\\{threeparttable\\}\\n\\\\begin\\{tabular\\}", perl = TRUE)
  tr <- texreg(model1, threeparttable = TRUE, longtable = TRUE, caption.above = TRUE)
  expect_match(tr, "\\usepackage\\{threeparttablex\\}\\n\\n")
  expect_match(tr, "\\\\hline\\n\\\\insertTableNotes\\\\\\\\\\n\\\\endlastfoot", perl = TRUE)
  expect_match(tr, "\\\\begin\\{TableNotes\\}\\[flushleft\\]\\n", perl = TRUE)
  expect_match(tr, "\\\\begin\\{ThreePartTable\\}\\n", perl = TRUE)
})

test_that("siunitx argument works in the texreg function", {
  tr <- texreg(list(model1, model1), single.row = TRUE, siunitx = TRUE, custom.gof.rows = list(abc = c("abcdefg hijklmn", "abc")))
  expect_match(tr, "\\& 3\\.31 \\\\; \\(0\\.06\\)\\^\\{\\*\\*\\*\\}  \\&", perl = TRUE)
  expect_match(tr, "l S\\[table-format=3\\.11\\]", perl = TRUE)
  expect_match(tr, "\\n\\\\usepackage\\{siunitx\\}\\n", perl = TRUE)
  expect_match(tr, "abc         \\& \\{abcdefg hijklmn\\}     \\& \\{abc\\}", perl = TRUE)
  expect_match(tr, "\\n\\\\sisetup\\{parse-numbers=false, table-text-alignment=right\\}\\n", perl = TRUE)
  tr <- texreg(list(model1, model1), siunitx = TRUE, custom.gof.rows = list(abc = c("abcdefg hijklmn", "abc")))
  expect_match(tr, "l S\\[table-format=3\\.5\\]", perl = TRUE)
  expect_match(tr, "\\{Model 1\\}", perl = TRUE)
  expect_warning(texreg(model1, siunitx = TRUE, dcolumn = TRUE), "The dcolumn and siunitx packages cannot be used at the same time. Switching off 'dcolumn'.")
  expect_warning(texreg(model1, siunitx = TRUE, bold = TRUE), "The siunitx package and the 'bold' argument cannot be used at the same time. Switching off 'siunitx'.")
  expect_warning(texreg(model1, dcolumn = TRUE, bold = TRUE), "The dcolumn package and the 'bold' argument cannot be used at the same time. Switching off 'dcolumn'.")
})

test_that("NA values in custom.coef.names are interpreted correctly", {
  mr <- matrixreg(model1, custom.coef.names = c("abc", NA))
  expect_match(mr[2, 1], "abc")
  expect_match(mr[4, 1], "Petal.Width")
  mr <- matrixreg(model1, custom.coef.names = c("abc", "abc"))
  expect_match(mr[2, 1], "abc")
  expect_match(mr[4, 1], "abc")
})

test_that("duplicate row labels in custom.coef.names are merged when feasible", {
  skip_if_not_installed("nlme")
  require("nlme")
  model.1 <- lme(distance ~ age, data = Orthodont, random = ~ 1)
  model.2 <- lme(distance ~ Sex, data = Orthodont, random = ~ 1)


  screenreg(list(model.1, model.2), custom.coef.names = c("Intercept", "sex", "sex"))

  expect_equivalent(dim(matrixreg(list(model.1, model.2))), c(12, 3))
  expect_equivalent(dim(matrixreg(list(model.1, model.2),
                                  custom.coef.names = c("Intercept", "sex", "sex"))),
                    c(10, 3))
})

test_that("LaTeX code is escaped in GOF labels", {
  skip_if_not_installed("lme4", minimum_version = "1.1.34")
  skip_if_not_installed("Matrix", minimum_version = "1.6.1")
  skip_on_ci()
  require("lme4")
  iris$species_variable <- iris$Species
  fit <- lmer(Sepal.Length ~ Sepal.Width + (1 | species_variable), data = iris)
  expect_match(texreg(fit),
               regexp = "Num. groups: species\\\\_variable",
               perl = TRUE)
  expect_match(screenreg(fit),
               regexp = "Num. groups: species_variable",
               perl = TRUE)
})

test_that("no stars, single.row, and dcolumn work together", {
  expect_match(texreg(model1,
                      single.row = TRUE,
                      dcolumn = TRUE,
                      stars = numeric()),
               regexp = "l D{\\)}{\\)}{9\\)0}}",
               perl = TRUE)
})

test_that("table and tabular arguments work in texreg function", {
  tr <- texreg(model1, table = FALSE)
  expect_match(tr, "^\\n\\\\begin\\{tabular")
  tr <- texreg(model1, table = FALSE, tabular = FALSE)
  expect_match(tr, "^\\n\\\\hline")
  expect_warning(texreg(model1, tabular = FALSE), "Setting 'table = FALSE'")
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
    expect_error(knitreg(list(model1, model1)), regexp = "knitreg requires the 'knitr' package to be installed")
  })
  with_mock(requireNamespace = function (package, ...) {
      ifelse(package == "rmarkdown", return(FALSE), return(TRUE))
    }, {
    expect_error(knitreg(list(model1, model1)), regexp = "knitreg requires the 'rmarkdown' package to be installed")
  })
  skip_if_not_installed("knitr", minimum_version = "1.22")
  require("knitr")
  skip_if_not_installed("rmarkdown", minimum_version = "1.12")
  require("rmarkdown")
  expect_match(knitreg(model1), "Petal.Width")

  # the following evaluates that knitreg chooses and outputs the expected format
  knitr::opts_knit$set(out.format = "markdown")

  with_mock("rmarkdown::all_output_formats" = function (input) {"html_document"},
            "knitr::current_input" = function () {NULL},
            expect_equivalent(knitreg(model1), htmlreg(model1, doctype = FALSE)))

  with_mock("rmarkdown::all_output_formats" = function (input) {"bookdown::html_document2"},
            "knitr::current_input" = function () {NULL},
            expect_equivalent(knitreg(model1), htmlreg(model1, doctype = FALSE)))

  with_mock("rmarkdown::all_output_formats" = function (input) {"pdf_document"},
            "knitr::current_input" = function () {NULL},
            expect_equivalent(knitreg(model1), texreg(model1, use.packages = FALSE)))

  with_mock("rmarkdown::all_output_formats" = function (input) {"bookdown::pdf_document2"},
            "knitr::current_input" = function () {NULL},
            expect_equivalent(knitreg(model1), texreg(model1, use.packages = FALSE)))

  with_mock("rmarkdown::all_output_formats" = function (input) {"bookdown::pdf_book"},
            "knitr::current_input" = function () {NULL},
            expect_equivalent(knitreg(model1), texreg(model1, use.packages = FALSE)))

  # formatting table to test word output in knitreg
  mr <- matrixreg(model1, output.type = "ascii", include.attributes = FALSE, trim = TRUE)
  colnames(mr) <- mr[1, ]
  mr <- mr[-1, ]

  with_mock("rmarkdown::all_output_formats" = function (input) {"word_document"},
            "knitr::current_input" = function () {NULL},
            expect_equivalent(knitreg(model1), knitr::kable(mr)))

  with_mock("rmarkdown::all_output_formats" = function (input) {"bookdown::word_document2"},
            "knitr::current_input" = function () {NULL},
            expect_equivalent(knitreg(model1), knitr::kable(mr)))
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

test_that("Single coef and no inference in custom extract method works", {
  # as per https://github.com/leifeld/texreg/issues/209
  my_obj <- lm(freeny)
  my_obj2 <- lm(y ~ lag.quarterly.revenue - 1, data = freeny)
  extract_my_lm <- function(model) {
    createTexreg(
      coef.names = names(coef(model)),
      coef = coef(model),
      se = numeric(0),
      pvalues = numeric(0),
    )
  }
  setMethod("extract", signature = className("lm", "stats"), definition = extract_my_lm)
  expect_no_error(screenreg(my_obj))
  expect_no_error(screenreg(my_obj2))
})