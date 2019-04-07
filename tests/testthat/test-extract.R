context("Extract methods")
library("texreg")

# dynlm (dynlm) ----
test_that("extract dynlm objects from the dynlm package", {
  skip_if_not_installed("dynlm")
  skip_if_not_installed("datasets")
  require("dynlm")
  set.seed(12345)
  data("UKDriverDeaths", package = "datasets")
  uk <- log10(UKDriverDeaths)
  dfm <- dynlm(uk ~ L(uk, 1) + L(uk, 12))
  tr <- extract(dfm, include.rmse = TRUE)
  expect_equivalent(tr@coef, c(0.18, 0.43, 0.51), tolerance = 1e-2)
  expect_equivalent(tr@se, c(0.158, 0.053, 0.056), tolerance = 1e-2)
  expect_equivalent(tr@pvalues, c(0.25, 0.00, 0.00), tolerance = 1e-2)
  expect_equivalent(tr@gof, c(0.68, 0.68, 180, 0.04), tolerance = 1e-2)
  expect_length(tr@gof.names, 4)
  expect_length(tr@coef, 3)
  expect_equivalent(which(tr@pvalues < 0.05), 2:3)
  expect_equivalent(which(tr@gof.decimal), c(1, 2, 4))
})

# feis (feisr) ----
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
  expect_equivalent(tr@coef, 0.056, tolerance = 1e-3)
  expect_equivalent(tr@se, 0.0234, tolerance = 1e-3)
  expect_equivalent(tr@pvalues, 0.0165, tolerance = 1e-3)
  expect_equivalent(tr@gof, c(0.002, 0.002, 3100, 268, 0.312), tolerance = 1e-3)
  expect_length(tr@gof.names, 5)
  tr2 <- extract(feis2.mod)
  expect_length(tr2@coef, 6)
  expect_length(which(tr2@pvalues < 0.05), 2)
  expect_length(which(tr2@gof.decimal), 3)
})

# ivreg (AER) ----
test_that("extract ivreg objects from the AER package", {
  skip_if_not_installed("AER")
  require("AER")
  set.seed(12345)
  data("CigarettesSW", package = "AER")
  CigarettesSW$rprice <- with(CigarettesSW, price / cpi)
  CigarettesSW$rincome <- with(CigarettesSW, income/population / cpi)
  CigarettesSW$tdiff <- with(CigarettesSW, (taxs - tax) / cpi)
  fm <- ivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff + I(tax/cpi),
              data = CigarettesSW,
              subset = year == "1995")
  tr1 <- extract(fm, vcov = sandwich, df = Inf, diagnostics = TRUE, include.rmse = TRUE)
  fm2 <- ivreg(log(packs) ~ log(rprice) | tdiff, data = CigarettesSW,
               subset = year == "1995")
  tr2 <- extract(fm2)
  expect_equivalent(tr1@coef, c(9.89, -1.28, 0.28), tolerance = 1e-2)
  expect_equivalent(tr1@se, c(0.93, 0.24, 0.25), tolerance = 1e-2)
  expect_equivalent(tr1@pvalues, c(0.00, 0.00, 0.25), tolerance = 1e-2)
  expect_equivalent(tr1@gof, c(0.43, 0.40, 48, 0.19), tolerance = 1e-2)
  expect_length(tr1@gof.names, 4)
  expect_length(tr2@coef, 2)
  expect_length(which(tr2@pvalues < 0.05), 2)
  expect_equivalent(which(tr2@gof.decimal), 1:2)
})

# lm (stats) ----
test_that("extract lm objects from the stats package", {
  set.seed(12345)
  ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
  trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
  group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
  weight <- c(ctl, trt)
  lm.D9 <- lm(weight ~ group)
  lm.D90 <- lm(weight ~ group - 1)
  tr <- extract(lm.D9)
  expect_equivalent(tr@coef, c(5.032, -0.371), tolerance = 1e-3)
  expect_equivalent(tr@se, c(0.22, 0.31), tolerance = 1e-2)
  expect_equivalent(tr@pvalues, c(0.00, 0.25), tolerance = 1e-2)
  expect_equivalent(tr@gof, c(0.07, 0.02, 20), tolerance = 1e-2)
  expect_length(tr@gof.names, 3)
  tr2 <- extract(lm.D90, include.rmse = TRUE)
  expect_length(tr2@coef, 2)
  expect_length(which(tr2@pvalues < 0.05), 2)
  expect_length(which(tr2@gof.decimal), 3)
})

# speedlm (speedglm) ----
test_that("extract speedlm objects from the speedglm package", {
  skip_if_not_installed("speedglm")
  require("speedglm")
  set.seed(12345)
  n <- 1000
  k <- 3
  y <- rnorm(n)
  x <- round(matrix(rnorm(n * k), n, k), digits = 3)
  colnames(x) <- c("s1", "s2", "s3")
  da <- data.frame(y, x)
  do1 <- da[1:300, ]
  do2 <- da[301:700, ]
  do3 <- da[701:1000, ]
  m1 <- speedlm(y ~ s1 + s2 + s3, data = do1)
  m1 <- update(m1, data = do2)
  m1 <- update(m1, data = do3)
  tr <- extract(m1, include.fstatistic = TRUE)
  expect_equivalent(tr@coef, c(0.05, 0.04, -0.01, -0.03), tolerance = 1e-2)
  expect_equivalent(tr@se, c(0.03, 0.03, 0.03, 0.03), tolerance = 1e-2)
  expect_equivalent(tr@pvalues, c(0.13, 0.22, 0.69, 0.39), tolerance = 1e-2)
  expect_equivalent(tr@gof, c(0, 0, 1000, 0.80), tolerance = 1e-2)
  expect_length(tr@gof.names, 4)
  expect_length(tr@coef, 4)
  expect_equivalent(which(tr@pvalues < 0.05), integer())
  expect_equivalent(which(tr@gof.decimal), c(1, 2, 4))
})