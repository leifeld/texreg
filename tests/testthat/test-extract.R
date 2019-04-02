context("Test of extract methods")
library("texreg")

test_that("extract feis objects from the feisr package", {
  skip_if_not_installed("feisr", minimum_version = "1.0.1")
  require("feisr")
  set.seed(12345)
  data("mwp", package = "feisr")
  feis1.mod <- feis(lnw ~ marry | exp, data = mwp, id = "id")
  feis2.mod <- feis(lnw ~ marry + enrol + as.factor(yeargr) | exp,
                    data = mwp,
                    id = "id")
  tr <- extract(feis1.mod)
  expect_equal(tr@coef, 0.056, tolerance = 1e-3)
  expect_equal(tr@se, 0.0234, tolerance = 1e-3)
  expect_equal(tr@pvalues, 0.0165, tolerance = 1e-3)
  expect_equivalent(tr@gof, c(0.002, 0.002, 3100, 268, 0.312), tolerance = 1e-3)
  expect_length(tr@gof.names, 5)
  tr2 <- extract(feis2.mod)
  expect_length(tr2@coef, 6)
  expect_length(which(tr2@pvalues < 0.05), 2)
  expect_length(which(tr2@gof.decimal), 3)
})