context("extract methods")
suppressPackageStartupMessages(library("texreg"))

# Arima (stats) ----
test_that("extract Arima objects from the stats package", {
  testthat::skip_on_cran()
  set.seed(12345)
  m <- arima(USAccDeaths,
             order = c(0, 1, 1),
             seasonal = list(order = c(0, 1, 1)))
  tr <- extract(m)
  expect_length(tr@coef.names, 2)
  expect_length(tr@coef, 2)
  expect_length(tr@se, 2)
  expect_length(tr@pvalues, 2)
  expect_length(tr@ci.low, 0)
  expect_length(tr@ci.up, 0)
  expect_length(tr@gof, 4)
  expect_length(tr@gof.names, 4)
  expect_length(tr@gof.decimal, 4)
  expect_equivalent(which(tr@gof.decimal), 1:3)
  expect_equivalent(which(tr@pvalues < 0.05), 1:2)
  expect_equivalent(dim(matrixreg(m)), c(9, 2))
})

# forecast_ARIMA (forecast) ----
test_that("extract forecast_ARIMA objects from the forecast package", {
  testthat::skip_on_cran()
  skip_if_not_installed("forecast")
  require("forecast")
  set.seed(12345)
  air.model <- Arima(window(AirPassengers, end = 1956 + 11 / 12),
                     order = c(0, 1, 1),
                     seasonal = list(order = c(0, 1, 1), period = 12),
                     lambda = 0)
  tr <- extract(air.model)
  expect_length(tr@coef.names, 2)
  expect_length(tr@coef, 2)
  expect_length(tr@se, 2)
  expect_length(tr@pvalues, 2)
  expect_length(tr@ci.low, 0)
  expect_length(tr@ci.up, 0)
  expect_length(tr@gof, 5)
  expect_length(tr@gof.names, 5)
  expect_length(tr@gof.decimal, 5)
  expect_equivalent(which(tr@gof.decimal), 1:4)
  expect_equivalent(which(tr@pvalues < 0.05), 1:2)
  expect_equivalent(dim(matrixreg(air.model)), c(10, 2))

  m1 <- arima(USAccDeaths,
              order = c(0, 1, 1),
              seasonal = list(order = c(0, 1, 1)))
  m2 <- Arima(USAccDeaths,
              order = c(0, 1, 1),
              seasonal = list(order = c(0, 1, 1)))
  expect_s3_class(m1, "Arima")
  expect_s3_class(m2, "Arima")
  expect_s3_class(m2, "forecast_ARIMA")
  m <- matrixreg(list(m1, m2))
  expect_equivalent(dim(m), c(10, 3))
  expect_equivalent(m[2:9, 2], m[2:9, 3])
  expect_equivalent(m[10, 1], "AICc")
})

# bergm (Bergm) ----
test_that("extract bergm objects from the Bergm package", {
  testthat::skip_on_cran()
  suppressWarnings(skip_if_not_installed("Bergm", minimum_version = "5.0.2"))
  require("Bergm")
  set.seed(12345)
  data(florentine)
  suppressWarnings(suppressMessages(
    p.flo <- bergm(flomarriage ~ edges + kstar(2),
                   burn.in    = 10,
                   aux.iters  = 30,
                   main.iters = 30,
                   gamma      = 1.2)))
  tr <- extract(p.flo)
  expect_length(tr@se, 0)
  expect_length(tr@pvalues, 0)
  expect_length(tr@ci.low, 2)
  expect_length(tr@ci.up, 2)
  expect_length(tr@gof, 0)
  expect_length(tr@coef, 2)
  expect_equivalent(dim(matrixreg(p.flo)), c(5, 2))
})

# bife (bife) ----
test_that("extract bife objects from the bife package", {
  testthat::skip_on_cran()
  skip_if_not_installed("bife", minimum_version = "0.7")
  require("bife")
  set.seed(12345)

  mod <- bife(LFP ~ I(AGE^2) + log(INCH) + KID1 + KID2 + KID3 + factor(TIME) | ID, psid)
  tr <- extract(mod)

  expect_length(tr@coef.names, 13)
  expect_length(tr@coef, 13)
  expect_length(tr@se, 13)
  expect_length(tr@pvalues, 13)
  expect_length(tr@ci.low, 0)
  expect_length(tr@ci.up, 0)
  expect_length(tr@gof, 3)
  expect_length(tr@gof.names, 3)
  expect_length(tr@gof.decimal, 3)
  expect_equivalent(which(tr@gof.decimal), 1:2)
  expect_equivalent(which(tr@pvalues < 0.05), c(1:4, 8:13))
  expect_equivalent(dim(matrixreg(mod)), c(30, 2))
})

## commented out because it takes long and causes segfault in combination with other tests
# # brmsfit (brms) ----
# test_that("extract brmsfit objects from the brms package", {
#   testthat::skip_on_cran()
#   skip_if_not_installed("brms", minimum_version = "2.8.8")
#   skip_if_not_installed("coda", minimum_version = "0.19.2")
#   require("brms")
#   require("coda")
#
#   # example 2 from brm help page; see ?brm
#   sink(nullfile())
#   suppressMessages(
#     fit2 <- brm(rating ~ period + carry + cs(treat),
#                 data = inhaler, family = sratio("logit"),
#                 prior = set_prior("normal(0,5)"), chains = 1))
#   sink()
#
#   suppressWarnings(tr <- extract(fit2))
#   expect_length(tr@gof.names, 4)
#   expect_length(tr@coef, 8)
#   expect_length(tr@se, 8)
#   expect_length(tr@pvalues, 0)
#   expect_length(tr@ci.low, 8)
#   expect_length(tr@ci.up, 8)
#   expect_equivalent(which(tr@gof.decimal), c(1, 3, 4))
#   suppressWarnings(expect_equivalent(dim(matrixreg(fit2)), c(21, 2)))
#
#   # example 1 from brm help page; see ?brm
#   bprior1 <- prior(student_t(5, 0, 10), class = b) + prior(cauchy(0, 2), class = sd)
#   sink(nullfile())
#   suppressMessages(
#     fit1 <- brm(count ~ zAge + zBase * Trt + (1|patient),
#                 data = epilepsy,
#                 family = poisson(),
#                 prior = bprior1))
#   sink()
#
#   expect_warning(suppressMessages(tr <- extract(fit1, use.HDI = TRUE, reloo = TRUE)))
#   expect_length(tr@gof.names, 5)
#   expect_length(tr@coef, 5)
#   expect_length(tr@se, 5)
#   expect_length(tr@pvalues, 0)
#   expect_length(tr@ci.low, 5)
#   expect_length(tr@ci.up, 5)
#   expect_equivalent(which(tr@gof.decimal), c(1:2, 4:5))
#   expect_equivalent(suppressWarnings(dim(matrixreg(fit1))), c(16, 2))
# })

# btergm (btergm) ----
test_that("extract btergm objects from the btergm package", {
  testthat::skip_on_cran()
  skip_if_not_installed("btergm", minimum_version = "1.10.10")
  set.seed(5)
  networks <- list()
  for (i in 1:10) {              # create 10 random networks with 10 actors
    mat <- matrix(rbinom(100, 1, .25), nrow = 10, ncol = 10)
    diag(mat) <- 0               # loops are excluded
    networks[[i]] <- mat         # add network to the list
  }

  covariates <- list()
  for (i in 1:10) {              # create 10 matrices as covariate
    mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
    covariates[[i]] <- mat       # add matrix to the list
  }

  suppressWarnings(fit <- btergm::btergm(networks ~ edges + istar(2) + edgecov(covariates), R = 100, verbose = FALSE))
  tr <- extract(fit)
  expect_length(tr@se, 0)
  expect_length(tr@pvalues, 0)
  expect_length(tr@ci.low, 3)
  expect_length(tr@ci.up, 3)
  expect_length(tr@gof, 1)
  expect_length(tr@coef, 3)
  expect_equivalent(dim(matrixreg(fit)), c(8, 2))
  expect_true(all(tr@ci.low < tr@coef))
  expect_true(all(tr@coef < tr@ci.up))
})

# clm (ordinal) ----
test_that("extract clm objects from the ordinal package", {
  testthat::skip_on_cran()
  skip_if_not_installed("ordinal", minimum_version = "2019.12.10")
  set.seed(12345)
  fit <- ordinal::clm(Species ~ Sepal.Length, data = iris)
  tr <- extract(fit)
  expect_length(tr@coef.names, 3)
  expect_length(tr@coef, 3)
  expect_length(tr@se, 3)
  expect_length(tr@pvalues, 3)
  expect_length(tr@ci.low, 0)
  expect_length(tr@ci.up, 0)
  expect_length(tr@gof, 4)
  expect_length(tr@gof.names, 4)
  expect_length(tr@gof.decimal, 4)
  expect_equivalent(which(tr@gof.decimal), 1:3)
  expect_equivalent(which(tr@pvalues < 0.05), 1:3)
  expect_equivalent(dim(matrixreg(fit)), c(11, 2))
})

# dynlm (dynlm) ----
test_that("extract dynlm objects from the dynlm package", {
  testthat::skip_on_cran()
  skip_if_not_installed("dynlm")
  skip_if_not_installed("datasets")
  require("dynlm")
  set.seed(12345)
  data("UKDriverDeaths", package = "datasets")
  uk <- log10(UKDriverDeaths)
  dfm <- dynlm(uk ~ L(uk, 1) + L(uk, 12))
  tr <- extract(dfm, include.rmse = TRUE)
  expect_length(tr@coef.names, 3)
  expect_length(tr@coef, 3)
  expect_length(tr@se, 3)
  expect_length(tr@pvalues, 3)
  expect_length(tr@ci.low, 0)
  expect_length(tr@ci.up, 0)
  expect_length(tr@gof, 4)
  expect_length(tr@gof.names, 4)
  expect_length(tr@gof.decimal, 4)
  expect_equivalent(which(tr@gof.decimal), c(1, 2, 4))
  expect_equivalent(which(tr@pvalues < 0.05), 2:3)
  expect_equivalent(dim(matrixreg(dfm)), c(10, 2))
})

# ergm (ergm) ----
test_that("extract ergm objects from the ergm package", {
  testthat::skip_on_cran()
  skip_if_not_installed("ergm", minimum_version = "4.1.2")
  require("ergm")

  set.seed(12345)
  data(florentine)
  suppressMessages(gest <- ergm(flomarriage ~ edges + absdiff("wealth")))
  tr1 <- extract(gest)
  expect_length(tr1@coef.names, 2)
  expect_length(tr1@coef, 2)
  expect_length(tr1@se, 2)
  expect_length(tr1@pvalues, 2)
  expect_length(tr1@ci.low, 0)
  expect_length(tr1@ci.up, 0)
  expect_length(tr1@gof, 3)
  expect_length(tr1@gof.names, 3)
  expect_length(tr1@gof.decimal, 3)
  expect_equivalent(which(tr1@gof.decimal), 1:3)
  expect_equivalent(dim(matrixreg(gest)), c(8, 2))

  data(molecule)
  molecule %v% "atomic type" <- c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3,
                                  3, 3, 3, 3, 3)
  suppressMessages(gest <- ergm(molecule ~ edges + kstar(2) + triangle +
                                  nodematch("atomic type")))
  tr2 <- extract(gest)
  expect_length(tr2@coef.names, 4)
  expect_length(tr2@coef, 4)
  expect_length(tr2@se, 4)
  expect_length(tr2@pvalues, 4)
  expect_length(tr2@ci.low, 0)
  expect_length(tr2@ci.up, 0)
  expect_length(tr2@gof, 3)
  expect_length(tr2@gof.names, 3)
  expect_length(tr2@gof.decimal, 3)
  expect_equivalent(which(tr2@gof.decimal), 1:3)
  expect_equivalent(dim(matrixreg(gest)), c(12, 2))
})

# feglm (alpaca) ----
test_that("extract feglm objects from the alpaca package", {
  testthat::skip_on_cran()
  skip_if_not_installed("alpaca", minimum_version = "0.3.2")
  require("alpaca")

  set.seed(12345)
  data <- simGLM(1000L, 20L, 1805L, model = "logit")
  mod <- feglm(y ~ x1 + x2 + x3 | i + t, data)

  tr <- extract(mod)
  expect_length(tr@coef.names, 3)
  expect_length(tr@coef, 3)
  expect_length(tr@se, 3)
  expect_length(tr@pvalues, 3)
  expect_length(tr@ci.low, 0)
  expect_length(tr@ci.up, 0)
  expect_length(tr@gof, 4)
  expect_length(tr@gof.names, 4)
  expect_length(tr@gof.decimal, 4)
  expect_equivalent(which(tr@gof.decimal), 1)
  expect_equivalent(which(tr@pvalues < 0.05), 1:3)
  expect_equivalent(dim(matrixreg(mod)), c(11, 2))
})

# feis (feisr) ----
test_that("extract feis objects from the feisr package", {
  testthat::skip_on_cran()
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

# felm (lfe) ----
test_that("extract felm objects from the lfe package", {
  testthat::skip_on_cran()
  skip_if_not_installed("lfe", minimum_version = "2.8.5")
  require("lfe")

  set.seed(12345)
  x <- rnorm(1000)
  x2 <- rnorm(length(x))
  id <- factor(sample(20, length(x), replace = TRUE))
  firm <- factor(sample(13, length(x),replace = TRUE))
  id.eff <- rnorm(nlevels(id))
  firm.eff <- rnorm(nlevels(firm))
  u <- rnorm(length(x))
  y <- x + 0.5 * x2 + id.eff[id] + firm.eff[firm] + u
  est <- felm(y ~ x + x2 | id + firm)

  tr <- extract(est)

  expect_equivalent(tr@coef, c(1.0188, 0.5182), tolerance = 1e-2)
  expect_equivalent(tr@se, c(0.032, 0.032), tolerance = 1e-2)
  expect_equivalent(tr@pvalues, c(0.00, 0.00), tolerance = 1e-2)
  expect_equivalent(tr@gof, c(1000, 0.7985, 0.575, 0.792, 0.560, 20, 13), tolerance = 1e-2)
  expect_length(tr@gof.names, 7)
  expect_length(tr@coef, 2)
  expect_equivalent(which(tr@pvalues < 0.05), 1:2)
  expect_equivalent(which(tr@gof.decimal), 2:5)

  # check exclusion of projected model statistics
  tr <- extract(est, include.proj.stats = FALSE)
  expect_length(tr@gof.names, 5)
  expect_false(any(grepl('proj model', tr@gof.names, fixed = TRUE)))

  # without fixed effects
  OLS1 <- felm(Sepal.Length ~ Sepal.Width |0|0|0, data = iris)
  tr1 <- extract(OLS1)
  expect_length(tr1@gof, 5)
})

# fixest (fixest) ----
test_that("extract fixest objects created with the fixest package", {
  testthat::skip_on_cran()
  skip_if_not_installed("fixest", minimum_version = "0.10.5")
  require("fixest")

  # test ordinary least squares with multiple fixed effects
  set.seed(12345)
  x <- rnorm(1000)
  data <- data.frame(
    x = x,
    x2 = rnorm(length(x)),
    id = factor(sample(20, length(x), replace = TRUE)),
    firm = factor(sample(13, length(x),replace = TRUE))
  )
  id.eff <- rnorm(nlevels(data$id))
  firm.eff <- rnorm(nlevels(data$firm))
  u <- rnorm(length(x))
  data$y <- with(data, x + 0.5 * x2 + id.eff[id] + firm.eff[firm] + u)
  est <- feols(y ~ x + x2 | id + firm, data = data)

  tr <- extract(est)

  expect_equivalent(tr@coef, c(1.0188, 0.5182), tolerance = 1e-2)
  # NOTE: standard errors differ from default produced by lfe (tested above)
  #       see https://cran.r-project.org/web/packages/fixest/vignettes/standard_errors.html
  expect_equivalent(tr@se, c(0.021, 0.032), tolerance = 1e-2)
  expect_equivalent(tr@pvalues, c(0.00, 0.00), tolerance = 1e-2)
  expect_equivalent(tr@gof, c(1000, 20, 13, 0.7985, 0.575, 0.792, 0.57), tolerance = 1e-2)
  expect_lte(length(tr@gof.names), 7)
  expect_gte(length(tr@gof.names), 5)
  expect_length(tr@coef, 2)
  expect_equivalent(which(tr@pvalues < 0.05), 1:2)

  # test generalized linear model
  data$y <- rpois(length(data$x), exp(data$x + data$x2 + id.eff[data$id]))
  est <- fepois(y ~ x + x2 | id, data = data)
  tr <- extract(est)

  expect_equivalent(tr@coef, c(1.00, 1.00), tolerance = 1e-2)
  expect_equivalent(tr@se, c(0.01, 0.02), tolerance = 1e-2)
  expect_equivalent(tr@pvalues, c(0.00, 0.00), tolerance = 1e-2)
  expect_equivalent(tr@gof, c(1000, 20, 955.4, -1479.6, 0.83), tolerance = 1e-2)
  expect_length(tr@gof.names, 5)
  expect_length(tr@coef, 2)
  expect_equivalent(which(!tr@gof.decimal), 1:2)
})

# gamlssZadj (gamlss.inf) ----
test_that("extract gamlssZadj objects from the gamlss.inf package", {
  testthat::skip_on_cran()
  skip_if_not_installed("gamlss.inf", minimum_version = "1.0.1")
  require("gamlss.inf")

  set.seed(12345)
  sink(nullfile())
  y0 <- rZAGA(1000, mu = .3, sigma = .4, nu = .15)
  g0 <- gamlss(y0 ~ 1, family = ZAGA)
  t0 <- gamlssZadj(y = y0, mu.formula = ~1, family = GA, trace = TRUE)
  sink()

  tr <- extract(t0)
  expect_length(tr@gof.names, 2)
  expect_length(tr@coef, 3)
  expect_length(tr@se, 3)
  expect_length(tr@pvalues, 3)
  expect_length(tr@ci.low, 0)
  expect_length(tr@ci.up, 0)
  expect_equivalent(which(tr@gof.decimal), 2)
  expect_equivalent(tr@coef.names, c("$\\mu$ (Intercept)",
                                     "$\\sigma$ (Intercept)",
                                     "$\\nu$ (Intercept)"))
})

# glm.cluster (miceadds) ----
test_that("extract glm.cluster objects from the miceadds package", {
  testthat::skip_on_cran()
  skip_if_not_installed("miceadds", minimum_version = "3.8.9")
  require("miceadds")

  data(data.ma01)
  dat <- data.ma01

  dat$highmath <- 1 * (dat$math > 600)
  mod2 <- miceadds::glm.cluster(data = dat,
                                formula = highmath ~ hisei + female,
                                cluster = "idschool",
                                family = "binomial")
  tr <- extract(mod2)

  expect_equivalent(tr@coef, c(-2.76, 0.03, -0.15), tolerance = 1e-2)
  expect_equivalent(tr@se, c(0.25, 0.00, 0.10), tolerance = 1e-2)
  expect_equivalent(tr@pvalues, c(0.00, 0.00, 0.13), tolerance = 1e-2)
  expect_equivalent(tr@gof, c(3108.095, 3126.432, -1551.047, 3102.095, 3336.000), tolerance = 1e-2)
  expect_length(tr@gof.names, 5)
  expect_length(tr@coef, 3)
  expect_equivalent(which(tr@pvalues < 0.05), 1:2)
  expect_equivalent(which(tr@gof.decimal), 1:4)
})

# glmerMod (lme4) ----
test_that("extract glmerMod objects from the lme4 package", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  skip_if_not_installed("lme4", minimum_version = "1.1.34")
  skip_if_not_installed("Matrix", minimum_version = "1.6.1")
  require("lme4")
  set.seed(12345)
  gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
               data = cbpp,
               family = binomial)
  expect_equivalent(class(gm1)[1], "glmerMod")
  tr <- extract(gm1, include.dic = TRUE, include.deviance = TRUE)
  expect_equivalent(tr@coef, c(-1.40, -0.99, -1.13, -1.58), tolerance = 1e-2)
  expect_equivalent(tr@se, c(0.23, 0.30, 0.32, 0.42), tolerance = 1e-2)
  expect_equivalent(tr@pvalues, c(0, 0, 0, 0), tolerance = 1e-2)
  expect_length(tr@gof.names, 8)
  expect_equivalent(which(tr@gof.decimal), c(1:5, 8))
  expect_length(which(grepl("Var", tr@gof.names)), 1)
  expect_length(which(grepl("Cov", tr@gof.names)), 0)
  tr_profile <- extract(gm1, method = "profile", nsim = 5)
  tr_boot <- suppressWarnings(extract(gm1, method = "boot", nsim = 5))
  tr_wald <- extract(gm1, method = "Wald")
  expect_length(tr_profile@se, 0)
  expect_length(tr_profile@ci.low, 4)
  expect_length(tr_profile@ci.up, 4)
  expect_length(tr_boot@se, 0)
  expect_length(tr_boot@ci.low, 4)
  expect_length(tr_boot@ci.up, 4)
  expect_length(tr_wald@se, 0)
  expect_length(tr_wald@ci.low, 4)
  expect_length(tr_wald@ci.up, 4)
})

# glmmTMB (glmmTMB) ----
test_that("extract glmmTMB objects from the glmmTMB package", {
  testthat::skip_on_cran()
  skip_if_not_installed("glmmTMB", minimum_version = "1.0.1")
  require("glmmTMB")

  set.seed(12345)
  m2 <- glmmTMB(count ~ spp + mined + (1|site),
                zi = ~ spp + mined,
                family = nbinom2, data = Salamanders)

  tr <- extract(m2)
  expect_length(tr@gof.names, 5)
  expect_length(tr@coef, 16)
  expect_length(tr@se, 16)
  expect_length(tr@pvalues, 16)
  expect_length(tr@ci.low, 0)
  expect_length(tr@ci.up, 0)
  expect_equivalent(which(tr@gof.decimal), c(1, 2, 5))

  tr <- extract(m2, beside = TRUE)
  expect_length(tr[[1]]@gof.names, 5)
  expect_length(tr[[1]]@coef, 8)
  expect_length(tr[[2]]@coef, 8)
  expect_length(tr[[1]]@se, 8)
  expect_length(tr[[2]]@se, 8)
  expect_length(tr[[1]]@pvalues, 8)
  expect_length(tr[[2]]@pvalues, 8)
  expect_length(tr, 2)
  expect_equivalent(which(tr[[2]]@gof.decimal), c(1, 2, 5))

  data("mtcars")
  cars <- glmmTMB(gear ~ mpg, data = mtcars)
  tr_cars <- extract(cars)
  expect_length(tr_cars@gof, 3)
  expect_equal(tr_cars@gof.decimal, c(TRUE, TRUE, FALSE))
  expect_equal(tr_cars@gof.names, c("AIC", "Log Likelihood", "Num. obs."))
  expect_length(tr_cars@coef, 2)
  expect_length(tr_cars@se, 2)
  expect_length(tr_cars@pvalues, 2)
})

# ivreg (AER) ----
test_that("extract ivreg objects from the AER package", {
  testthat::skip_on_cran()
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

# lm.cluster (miceadds) ----
test_that("extract lm.cluster objects from the miceadds package", {
  testthat::skip_on_cran()
  skip_if_not_installed("miceadds", minimum_version = "3.8.9")
  require("miceadds")

  data(data.ma01)
  dat <- data.ma01

  mod1 <- miceadds::lm.cluster(data = dat,
                               formula = read ~ hisei + female,
                               cluster = "idschool")
  tr <- extract(mod1)

  expect_equivalent(tr@coef, c(418.80, 1.54, 35.70), tolerance = 1e-2)
  expect_equivalent(tr@se, c(6.45, 0.11, 3.81), tolerance = 1e-2)
  expect_equivalent(tr@pvalues, c(0.00, 0.00, 0.00), tolerance = 1e-2)
  expect_equivalent(tr@gof, c(0.15, 0.15, 3180), tolerance = 1e-2)
  expect_length(tr@gof.names, 3)
  expect_length(tr@coef, 3)
  expect_equivalent(which(tr@pvalues < 0.05), 1:3)
  expect_equivalent(which(tr@gof.decimal), 1:2)
})

# lmerMod (lme4) ----
test_that("extract lmerMod objects from the lme4 package", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  skip_if_not_installed("lme4", minimum_version = "1.1.34")
  skip_if_not_installed("Matrix", minimum_version = "1.6.1")
  require("lme4")
  set.seed(12345)
  fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  fm1_ML <- update(fm1, REML = FALSE)
  fm2 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)
  tr1 <- extract(fm1, include.dic = TRUE, include.deviance = TRUE)
  tr1_ML <- extract(fm1_ML, include.dic = TRUE, include.deviance = TRUE)
  tr2_profile <- extract(fm2, method = "profile", nsim = 5)
  tr2_boot <- suppressWarnings(extract(fm2, method = "boot", nsim = 5))
  tr2_wald <- extract(fm2, method = "Wald")
  expect_equivalent(class(fm1)[1], "lmerMod")
  expect_equivalent(tr1@coef, c(251.41, 10.47), tolerance = 1e-2)
  expect_equivalent(tr1@coef, tr1_ML@coef, tolerance = 1e-2)
  expect_equivalent(tr1@se, c(6.82, 1.55), tolerance = 1e-2)
  expect_equivalent(tr1@pvalues, c(0, 0), tolerance = 1e-2)
  expect_equivalent(tr1@gof, c(1755.63, 1774.79, 1760.25, 1751.94, -871.81, 180, 18, 611.90, 35.08, 9.61, 654.94), tolerance = 1e-2)
  expect_length(tr1@gof.names, 11)
  expect_equivalent(which(tr1@gof.decimal), c(1:5, 8:11))
  expect_equivalent(tr1@coef, tr1_ML@coef)
  expect_length(tr1_ML@gof, 11)
  expect_length(tr2_profile@gof, 8)
  expect_equivalent(tr1@coef, tr2_profile@coef, tolerance = 1e-2)
  expect_equivalent(tr1@coef, tr2_boot@coef, tolerance = 1e-2)
  expect_equivalent(tr1@coef, tr2_wald@coef, tolerance = 1e-2)
  expect_length(which(grepl("Var", tr1@gof.names)), 3)
  expect_length(which(grepl("Var", tr2_wald@gof.names)), 3)
  expect_length(which(grepl("Cov", tr1@gof.names)), 1)
  expect_length(which(grepl("Cov", tr2_wald@gof.names)), 0)
})

# maxLik (maxLik) ----
test_that("extract maxLik objects from the maxLik package", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("maxLik", minimum_version = "1.4.8")
  require("maxLik")
  set.seed(12345)

  # example 1 from help page
  t <- rexp(100, 2)
  loglik <- function(theta) log(theta) - theta * t
  gradlik <- function(theta) 1 / theta - t
  hesslik <- function(theta) -100 / theta^2
  sink(nullfile())
  a <- maxLik(loglik, start = 1, control = list(printLevel = 2))
  sink()

  tr1 <- extract(a)
  expect_length(tr1@coef.names, 1)
  expect_length(tr1@coef, 1)
  expect_length(tr1@se, 1)
  expect_length(tr1@pvalues, 1)
  expect_length(tr1@ci.low, 0)
  expect_length(tr1@ci.up, 0)
  expect_true(!any(is.na(tr1@coef)))
  expect_length(tr1@gof, 2)
  expect_length(tr1@gof.names, 2)
  expect_length(tr1@gof.decimal, 2)
  expect_equivalent(which(tr1@gof.decimal), 1:2)

  # example 2 from help page
  b <- maxLik(loglik, gradlik, hesslik, start = 1,
              control = list(tol = -1, reltol = 1e-12, gradtol = 1e-12))

  tr2 <- extract(b)
  expect_length(tr2@coef.names, 1)
  expect_length(tr2@coef, 1)
  expect_length(tr2@se, 1)
  expect_length(tr2@pvalues, 1)
  expect_length(tr2@ci.low, 0)
  expect_length(tr2@ci.up, 0)
  expect_true(!any(is.na(tr2@coef)))
  expect_length(tr2@gof, 2)
  expect_length(tr2@gof.names, 2)
  expect_length(tr2@gof.decimal, 2)
  expect_equivalent(which(tr2@gof.decimal), 1:2)

  # example 3 from help page
  loglik <- function(param) {
    mu <- param[1]
    sigma <- param[2]
    ll <- -0.5 * N * log(2 * pi) - N * log(sigma) - sum(0.5 * (x - mu)^2 / sigma^2)
    ll
  }
  x <- rnorm(100, 1, 2)
  N <- length(x)
  res <- maxLik(loglik, start = c(0, 1))

  tr3 <- extract(res)
  expect_length(tr3@coef.names, 2)
  expect_length(tr3@coef, 2)
  expect_length(tr3@se, 2)
  expect_length(tr3@pvalues, 2)
  expect_length(tr3@ci.low, 0)
  expect_length(tr3@ci.up, 0)
  expect_true(!any(is.na(tr3@coef)))
  expect_length(tr3@gof, 2)
  expect_length(tr3@gof.names, 2)
  expect_length(tr3@gof.decimal, 2)
  expect_equivalent(which(tr3@gof.decimal), 1:2)

  # example 4 from help page
  resFix <- maxLik(loglik, start = c(mu = 0, sigma = 1), fixed = "sigma")

  tr4 <- extract(resFix)
  expect_length(tr3@coef.names, 2)
  expect_length(tr3@coef, 2)
  expect_length(tr3@se, 2)
  expect_length(tr3@pvalues, 2)
  expect_length(tr3@ci.low, 0)
  expect_length(tr3@ci.up, 0)
  expect_true(!any(is.na(tr3@coef)))
  expect_length(tr3@gof, 2)
  expect_length(tr3@gof.names, 2)
  expect_length(tr3@gof.decimal, 2)
  expect_equivalent(which(tr3@gof.decimal), 1:2)
})

# mlogit (mlogit) ----
test_that("extract mlogit objects from the mlogit package", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("mlogit", minimum_version = "1.1.0")
  require("mlogit")
  set.seed(12345)
  data("Fishing", package = "mlogit")
  Fish <- dfidx(Fishing, varying = 2:9, shape = "wide", choice = "mode")
  m <- mlogit(mode ~ price + catch | income, data = Fish)
  tr1 <- extract(m)

  expect_equivalent(sum(abs(tr1@coef)), 3.382753, tolerance = 1e-2)
  expect_equivalent(sum(tr1@se), 0.7789933, tolerance = 1e-2)
  expect_equivalent(sum(tr1@pvalues), 0.6136796, tolerance = 1e-2)
  expect_equivalent(sum(tr1@gof), 2417.138, tolerance = 1e-2)
  expect_length(tr1@coef, 8)
  expect_length(tr1@gof, 4)
  expect_equivalent(which(tr1@gof.decimal), 1:2)
  expect_equivalent(tr1@gof[4], 4)
  expect_equal(dim(matrixreg(tr1)), c(21, 2))
  expect_warning(extract(m, beside = TRUE), "choice-specific covariates")
})

# multinom (nnet) ----
test_that("extract multinom objects from the nnet package", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("nnet", minimum_version = "7.3.12")
  require("nnet")

  # example from https://thomasleeper.com/Rcourse/Tutorials/nominalglm.html
  set.seed(100)
  y <- sort(sample(1:3, 600, TRUE))
  x <- numeric(length = 600)
  x[1:200] <- -1 * x[1:200] + rnorm(200, 4, 2)
  x[201:400] <- 1 * x[201:400] + rnorm(200)
  x[401:600] <- 2 * x[401:600] + rnorm(200, 2, 2)

  sink(nullfile())
  m1 <- multinom(y ~ x)
  sink()
  tr2 <- extract(m1, beside = FALSE)
  tr3 <- extract(m1, beside = TRUE)

  expect_equivalent(sum(abs(tr2@coef)), 6.845567, tolerance = 1e-2)
  expect_equivalent(sum(tr2@se), 0.6671602, tolerance = 1e-2)
  expect_equivalent(sum(tr2@pvalues), 1.677308e-16, tolerance = 1e-2)
  expect_equivalent(sum(tr2@gof), 2852.451, tolerance = 1e-2)
  expect_length(tr2@coef, 4)
  expect_length(tr2@gof, 6)
  expect_equivalent(which(tr2@gof.decimal), 1:4)
  expect_equivalent(tr2@gof[6], 3)
  expect_equal(dim(matrixreg(tr2)), c(15, 2))
  expect_length(tr3, 2)
  expect_length(tr3[[1]]@coef, 2)
  expect_length(tr3[[2]]@coef, 2)
})

# nlmerMod (lme4) ----
test_that("extract nlmerMod objects from the lme4 package", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  skip_if_not_installed("lme4", minimum_version = "1.1.34")
  skip_if_not_installed("Matrix", minimum_version = "1.6.1")
  require("lme4")
  set.seed(12345)
  startvec <- c(Asym = 200, xmid = 725, scal = 350)
  nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
               Orange,
               start = startvec)
  expect_equivalent(class(nm1)[1], "nlmerMod")
  expect_warning(extract(nm1, include.dic = TRUE, include.deviance = TRUE),
                 "falling back to var-cov estimated from RX")
  tr <- suppressWarnings(extract(nm1, include.dic = TRUE, include.deviance = TRUE))
  expect_equivalent(tr@coef, c(192.05, 727.90, 348.07), tolerance = 1e-2)
  expect_equivalent(tr@se, c(15.58, 34.44, 26.31), tolerance = 1e-2)
  expect_equivalent(tr@pvalues, c(0, 0, 0), tolerance = 1e-2)
  expect_length(tr@gof.names, 9)
  expect_equivalent(which(tr@gof.decimal), c(1:5, 8, 9))
  expect_length(which(grepl("Var", tr@gof.names)), 2)
  expect_length(which(grepl("Cov", tr@gof.names)), 0)
  tr_wald <- suppressWarnings(extract(nm1, method = "Wald"))
  expect_length(tr_wald@se, 0)
  expect_length(tr_wald@ci.low, 3)
  expect_length(tr_wald@ci.up, 3)
})

# pcce (plm) ----
test_that("extract pcce objects from the plm package", {
  testthat::skip_on_cran()
  skip_if_not_installed("plm", minimum_version = "2.4.1")
  require("plm")
  set.seed(12345)
  data("Produc", package = "plm")

  ccepmod <- pcce(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, data = Produc, model="p")
  tr <- extract(ccepmod)
  expect_length(tr@coef.names, 4)
  expect_length(tr@coef, 4)
  expect_length(tr@se, 4)
  expect_length(tr@pvalues, 4)
  expect_length(tr@ci.low, 0)
  expect_length(tr@ci.up, 0)
  expect_length(tr@gof, 4)
  expect_length(tr@gof.names, 4)
  expect_length(tr@gof.decimal, 4)
  expect_equivalent(which(tr@gof.decimal), 1:3)

  ccemgmod <- pcce(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, data = Produc, model="mg")
  tr2 <- extract(ccemgmod)
  expect_length(tr2@coef.names, 4)
  expect_length(tr2@coef, 4)
  expect_length(tr2@se, 4)
  expect_length(tr2@pvalues, 4)
  expect_length(tr2@ci.low, 0)
  expect_length(tr2@ci.up, 0)
  expect_length(tr2@gof, 4)
  expect_length(tr2@gof.names, 4)
  expect_length(tr2@gof.decimal, 4)
  expect_equivalent(which(tr2@gof.decimal), 1:3)
})

# prais (prais) ----
test_that("extract prais objects from the prais package", {
  testthat::skip_on_cran()
  skip_if_not_installed("prais", minimum_version = "1.1.3")
  set.seed(12345)
  data("barium", package = "prais")

  output <- capture.output({
    pwmod <- prais::prais_winsten(lchnimp ~ lchempi + lgas + lrtwex + befile6 + affile6 + afdec6, data = barium, index = "t")
  }, file = NULL)

  tr <- extract(pwmod)
  expect_length(tr@coef.names, 8)
  expect_length(tr@coef, 8)
  expect_length(tr@se, 8)
  expect_length(tr@pvalues, 8)
  expect_length(tr@ci.low, 0)
  expect_length(tr@ci.up, 0)
  expect_length(tr@gof, 3)
  expect_length(tr@gof.names, 3)
  expect_length(tr@gof.decimal, 3)
  expect_equivalent(which(tr@gof.decimal), 1:2)
})

# remstimate (remstimate) ----
test_that("extract remstimate objects from the remstimate package", {
  testthat::skip_on_cran()
  skip_if_not_installed("remstimate", minimum_version = "2.3.11")
  skip_if_not_installed("remify", minimum_version = "3.2.6")
  skip_if_not_installed("remstats", minimum_version = "3.2.2")

  data(tie_data, package = "remstimate")
  tie_reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")
  tie_model <- ~ 1 +
    remstats::indegreeSender() +
    remstats::inertia() +
    remstats::reciprocity()
  tie_reh_stats <- remstats::remstats(reh = tie_reh, tie_effects = tie_model)
  rem1 <- remstimate::remstimate(reh = tie_reh,
                                 stats = tie_reh_stats,
                                 method = "MLE",
                                 ncores = 1)
  rem2 <- remstimate::remstimate(reh = tie_reh,
                                 stats = tie_reh_stats,
                                 method = "HMC",
                                 ncores = 1,
                                 L = 5L)
  rem3 <- remstimate::remstimate(reh = tie_reh,
                                 stats = tie_reh_stats,
                                 method = "GDADAMAX",
                                 ncores = 1)
  rem4 <- remstimate::remstimate(reh = tie_reh,
                                 stats = tie_reh_stats,
                                 method = "BSIR",
                                 ncores = 1)
  mr1 <- matrixreg(list(rem1, rem2, rem3, rem4))
  expect_true("matrix" %in% class(mr1))
  expect_equal(nrow(mr1), 14)
  expect_equal(ncol(mr1), 5)

  actor_reh <- remify::remify(edgelist = tie_data$edgelist, model = "actor")
  sender_model <- ~ 1 + remstats::outdegreeSender()
  receiver_model <- ~ 1 + remstats::otp()
  actor_reh_stats <- remstats::remstats(reh = actor_reh, sender_effects = sender_model, receiver_effects = receiver_model)
  rem5 <- remstimate::remstimate(reh = actor_reh,
                                 stats = actor_reh_stats,
                                 method = "MLE",
                                 ncores = 1)
  rem6 <- remstimate::remstimate(reh = actor_reh,
                                 stats = actor_reh_stats,
                                 method = "HMC",
                                 ncores = 1,
                                 L = 5L)
  rem7 <- remstimate::remstimate(reh = actor_reh,
                                 stats = actor_reh_stats,
                                 method = "GDADAMAX",
                                 ncores = 1)
  tr5 <- extract(rem5)
  expect_length(tr5, 2)
  expect_length(tr5[[1]]@coef.names, 2)
  expect_length(tr5[[1]]@gof, 5)
  expect_equal(tr5[[1]]@model.name, "sender_model")
  expect_length(tr5[[2]]@coef.names, 1)
  expect_length(tr5[[2]]@gof, 5)
  expect_equal(tr5[[2]]@model.name, "receiver_model")
  mr2 <- matrixreg(list(rem5, rem6, rem7))
  expect_true("matrix" %in% class(mr2))
  expect_equal(nrow(mr2), 12)
  expect_equal(ncol(mr2), 7)
})

# Sarlm (spatialreg) ----
test_that("extract Sarlm objects from the spatialreg package", {
  testthat::skip_on_cran()
  skip_if_not_installed("spatialreg", minimum_version = "1.2.1")
  require("spatialreg")
  set.seed(12345)

  # first example from ?lagsarlm
  data(oldcol, package = "spdep")
  listw <- spdep::nb2listw(COL.nb, style = "W")
  ev <- spatialreg::eigenw(listw)
  W <- as(listw, "CsparseMatrix")
  trMatc <- spatialreg::trW(W, type = "mult")

  output <- capture.output({
    COL.lag.eig <- spatialreg::lagsarlm(CRIME ~ INC + HOVAL,
                                        data = COL.OLD,
                                        listw = listw,
                                        method = "eigen",
                                        quiet = FALSE,
                                        control = list(pre_eig = ev,
                                                       OrdVsign = 1))
  }, file = NULL)

  tr <- extract(COL.lag.eig)
  expect_length(tr@coef.names, 4)
  expect_length(tr@coef, 4)
  expect_length(tr@se, 4)
  expect_length(tr@pvalues, 4)
  expect_length(tr@ci.low, 0)
  expect_length(tr@ci.up, 0)
  expect_length(tr@gof, 7)
  expect_length(tr@gof.names, 7)
  expect_length(tr@gof.decimal, 7)
  expect_equivalent(which(tr@gof.decimal), 3:7)

  # example from ?predict.Sarlm
  lw <- spdep::nb2listw(COL.nb)
  COL.lag.eig2 <- COL.mix.eig <- lagsarlm(CRIME ~ INC + HOVAL,
                                          data = COL.OLD,
                                          lw,
                                          type = "mixed")
  tr2 <- extract(COL.lag.eig2)
  expect_length(tr2@coef.names, 6)
  expect_length(tr2@coef, 6)
  expect_length(tr2@se, 6)
  expect_length(tr2@pvalues, 6)
  expect_length(tr2@ci.low, 0)
  expect_length(tr2@ci.up, 0)
  expect_length(tr2@gof, 7)
  expect_length(tr2@gof.names, 7)
  expect_length(tr2@gof.decimal, 7)
  expect_equivalent(which(tr2@gof.decimal), 3:7)
})

# speedglm (speedglm) ----
test_that("extract speedglm objects from the speedglm package", {
  testthat::skip_on_cran()
  skip_if_not_installed("speedglm", minimum_version = "0.3.2")
  require("speedglm")
  set.seed(12345)
  n <- 50000
  k <- 80
  y <- rgamma(n, 1.5, 1)
  x <-round( matrix(rnorm(n * k), n, k), digits = 3)
  colnames(x) <-paste("s", 1:k, sep = "")
  da <- data.frame(y, x)
  fo <- as.formula(paste("y ~", paste(paste("s", 1:k, sep = ""), collapse = " + ")))
  m3 <- speedglm(fo, data = da, family = Gamma(log))
  tr <- extract(m3)
  expect_length(tr@gof.names, 5)
  expect_length(tr@coef, 81)
  expect_equivalent(tr@gof.names, c("AIC", "BIC", "Log Likelihood", "Deviance", "Num. obs."))
  expect_equivalent(which(tr@pvalues < 0.05), c(1, 4, 5, 17, 20, 21, 43, 65, 68, 73, 80))
  expect_equivalent(which(tr@gof.decimal), 1:4)
})

# speedlm (speedglm) ----
test_that("extract speedlm objects from the speedglm package", {
  testthat::skip_on_cran()
  skip_if_not_installed("speedglm", minimum_version = "0.3.2")
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

# truncreg (truncreg) ----
test_that("extract truncreg objects from the truncreg package", {
  testthat::skip_on_cran()
  skip_if_not_installed("truncreg", minimum_version = "0.2.5")
  require("truncreg")

  set.seed(12345)
  x <- rnorm(100, mean = 1)
  y <- rnorm(100, mean = 1.3)
  dta <- data.frame(x, y)
  dta <- dta[y < quantile(y, 0.8), ]
  model <- truncreg(y ~ x, data = dta, point = max(dta$y), direction = "right")
  tr <- extract(model)

  expect_equivalent(tr@coef, c(1.24, 0.05, 0.96), tolerance = 1e-2)
  expect_equivalent(tr@se, c(0.25, 0.12, 0.14), tolerance = 1e-2)
  expect_equivalent(tr@pvalues, c(0, 0.67, 0), tolerance = 1e-2)
  expect_equivalent(tr@gof, c(80, -81.69, 169.38, 176.53), tolerance = 1e-2)
  expect_length(tr@gof.names, 4)
  expect_length(tr@coef, 3)
  expect_equivalent(which(tr@pvalues < 0.05), c(1, 3))
  expect_equivalent(which(tr@gof.decimal), 2:4)
})

# weibreg (eha) ----
test_that("extract weibreg objects from the eha package", {
  testthat::skip_on_cran()
  skip_if_not_installed("eha", minimum_version = "2.9.0")
  require("eha")

  set.seed(12345)
  # stratified model example from weibreg help page
  dat <- data.frame(time = c(4, 3, 1, 1, 2, 2, 3),
                    status = c(1, 1, 1, 0, 1, 1, 0),
                    x = c(0, 2, 1, 1, 1, 0, 0),
                    sex = c(0, 0, 0, 0, 1, 1, 1))
  model <- eha::weibreg(Surv(time, status) ~ x + strata(sex), data = dat)
  tr <- extract(model)

  expect_length(tr@coef, 5)
  expect_equivalent(class(tr@coef), "numeric")
  expect_length(tr@se, 5)
  expect_equivalent(class(tr@se), "numeric")
  expect_length(tr@pvalues, 5)
  expect_equivalent(class(tr@pvalues), "numeric")
  expect_length(tr@coef.names, 5)
  expect_length(tr@ci.low, 0)
  expect_length(tr@ci.up, 0)
  expect_length(tr@gof, 6)
  expect_length(tr@gof.names, 6)
  expect_length(tr@gof.decimal, 6)
  expect_equivalent(tr@gof[5], 5)
  expect_equivalent(which(tr@pvalues < 0.05), 2:5)
  expect_equivalent(which(tr@gof.decimal), 1:3)
})

# wls (metaSEM) ----
test_that("extract wls objects from the metaSEM package", {
  testthat::skip_on_cran()
  skip_if_not_installed("metaSEM", minimum_version = "1.2.5.1")
  require("metaSEM")
  set.seed(12345)

  # example 1 from wls help page: analysis of correlation structure
  R1.labels <- c("a1", "a2", "a3", "a4")
  R1 <- matrix(c(1.00, 0.22, 0.24, 0.18,
                 0.22, 1.00, 0.30, 0.22,
                 0.24, 0.30, 1.00, 0.24,
                 0.18, 0.22, 0.24, 1.00), ncol = 4, nrow = 4,
               dimnames = list(R1.labels, R1.labels))
  n <- 1000
  acovR1 <- metaSEM::asyCov(R1, n)
  model1 <- "f =~ a1 + a2 + a3 + a4"
  RAM1 <- metaSEM::lavaan2RAM(model1, obs.variables = R1.labels)
  wls.fit1a <- metaSEM::wls(Cov = R1, aCov = acovR1, n = n, RAM = RAM1,
                            cor.analysis = TRUE, intervals = "LB")
  tr1 <- extract(wls.fit1a)
  expect_length(tr1@coef.names, 4)
  expect_length(tr1@coef, 4)
  expect_length(tr1@se, 0)
  expect_length(tr1@pvalues, 0)
  expect_length(tr1@ci.low, 4)
  expect_length(tr1@ci.up, 4)
  expect_true(!any(is.na(tr1@coef)))
  expect_length(tr1@gof, 11)
  expect_length(tr1@gof.names, 11)
  expect_length(tr1@gof.decimal, 11)
  expect_equivalent(tr1@gof[8], 0.23893943, tolerance = 1e-2)
  expect_equivalent(which(tr1@gof.decimal), c(1, 3, 4, 5, 6, 7, 8, 10, 11))

  # example 2 from wls help page: multiple regression
  R2.labels <- c("y", "x1", "x2")
  R2 <- matrix(c(1.00, 0.22, 0.24,
                 0.22, 1.00, 0.30,
                 0.24, 0.30, 1.00), ncol = 3, nrow = 3,
               dimnames = list(R2.labels, R2.labels))
  acovR2 <- metaSEM::asyCov(R2, n)
  model2 <- "y ~ x1 + x2
             ## Variances of x1 and x2 are 1
             x1 ~~ 1*x1
             x2 ~~ 1*x2
             ## x1 and x2 are correlated
             x1 ~~ x2"
  RAM2 <- metaSEM::lavaan2RAM(model2, obs.variables = R2.labels)
  wls.fit2a <- metaSEM::wls(Cov = R2, aCov = acovR2, n = n, RAM = RAM2,
                            cor.analysis = TRUE, intervals = "LB")
  tr2 <- extract(wls.fit2a)
  expect_length(tr2@coef.names, 3)
  expect_length(tr2@coef, 3)
  expect_length(tr2@se, 0)
  expect_length(tr2@pvalues, 0)
  expect_length(tr2@ci.low, 3)
  expect_length(tr2@ci.up, 3)
  expect_true(!any(is.na(tr2@coef)))
  expect_length(tr2@gof, 11)
  expect_length(tr2@gof.names, 11)
  expect_length(tr2@gof.decimal, 11)
  expect_equivalent(tr2@gof[8], 0.0738, tolerance = 1e-2)
  expect_equivalent(which(tr2@gof.decimal), c(1, 3, 4, 5, 6, 7, 8, 10, 11))

  # example 3 from wls help page
  R3.labels <- c("a1", "a2", "a3", "a4")
  R3 <- matrix(c(1.50, 0.22, 0.24, 0.18,
                 0.22, 1.60, 0.30, 0.22,
                 0.24, 0.30, 1.80, 0.24,
                 0.18, 0.22, 0.24, 1.30), ncol = 4, nrow = 4,
               dimnames = list(R3.labels, R3.labels))
  n <- 1000
  acovS3 <- metaSEM::asyCov(R3, n, cor.analysis = FALSE)
  model3 <- "f =~ a1 + a2 + a3 + a4"
  RAM3 <- metaSEM::lavaan2RAM(model3, obs.variables = R3.labels)
  wls.fit3a <- metaSEM::wls(Cov = R3, aCov = acovS3, n = n, RAM = RAM3,
                            cor.analysis = FALSE)
  tr3 <- extract(wls.fit3a)
  expect_length(tr3@coef.names, 8)
  expect_length(tr3@coef, 8)
  expect_length(tr3@se, 8)
  expect_length(tr3@pvalues, 8)
  expect_length(tr3@ci.low, 0)
  expect_length(tr3@ci.up, 0)
  expect_true(!any(is.na(tr3@coef)))
  expect_length(tr3@gof, 10)
  expect_length(tr3@gof.names, 10)
  expect_length(tr3@gof.decimal, 10)
  expect_equivalent(which(tr3@gof.decimal), c(1, 3, 4, 5, 6, 7, 9, 10))
  expect_true(all(tr3@pvalues < 0.05))
})

# logitr (logitr) ----
test_that("extract logitr objects from the logitr package", {
  testthat::skip_on_cran()
  skip_if_not_installed("logitr", minimum_version = "0.8.0")
  require("logitr")
  set.seed(12345)

  mnl_pref <- logitr(
    data    = yogurt,
    outcome = "choice",
    obsID   = "obsID",
    pars    = c("price", "feat", "brand")
  )
  tr <- extract(mnl_pref)

  expect_equivalent(tr@coef, c(-0.37,  0.49, -3.72, -0.64, 0.73), tolerance = 1e-2)
  expect_equivalent(tr@se, c(0.02, 0.12, 0.15, 0.05, 0.08), tolerance = 1e-2)
  expect_equivalent(tr@pvalues, c(0, 0, 0, 0, 0), tolerance = 1e-2)
  expect_equivalent(tr@gof, c(2412.00, -2656.89, 5323.78, 5352.72), tolerance = 1e-2)
  expect_equivalent(which(tr@gof.decimal), c(2, 3, 4))
  expect_equivalent(which(tr@pvalues < 0.05), seq(1, 5))
  expect_length(tr@coef.names, 5)
  expect_length(tr@coef, 5)
  expect_length(tr@se, 5)
  expect_length(tr@pvalues, 5)
  expect_length(tr@ci.low, 0)
  expect_length(tr@ci.up, 0)
  expect_length(tr@gof.names, 4)
  expect_length(tr@gof, 4)
  expect_length(tr@gof.decimal, 4)
  expect_equivalent(dim(matrixreg(mnl_pref)), c(15, 2))
})
