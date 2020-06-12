context("extract methods")
suppressPackageStartupMessages(library("texreg"))


# Arima (stats) ----
test_that("extract Arima objects from the stats package", {
  set.seed(12345)
  m <- arima(USAccDeaths,
             order = c(0, 1, 1),
             seasonal = list(order = c(0, 1, 1)))
  tr <- extract(m)
  expect_equivalent(tr@coef, c(-0.43, -0.55), tolerance = 1e-2)
  expect_equivalent(tr@se, c(0.12, 0.178), tolerance = 1e-2)
  expect_equivalent(tr@pvalues, c(0.00, 0.00), tolerance = 1e-2)
  expect_equivalent(tr@gof, c(856.88, 863.1126, -425.44, 59), tolerance = 1e-2)
  expect_length(tr@gof.names, 4)
  expect_length(tr@coef, 2)
  expect_equivalent(which(tr@pvalues < 0.05), 1:2)
  expect_equivalent(which(tr@gof.decimal), 1:3)
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
  expect_equivalent(tr@coef, c(-0.39, -0.61), tolerance = 1e-2)
  expect_equivalent(tr@se, c(0.117, 0.1075829), tolerance = 1e-2)
  expect_equivalent(tr@pvalues, c(0.00, 0.00), tolerance = 1e-2)
  expect_equivalent(tr@gof, c(-291.5260, -291.2222, -284.2695, 148.7630, 83), tolerance = 1e-2)
  expect_length(tr@gof.names, 5)
  expect_length(tr@coef, 2)
  expect_equivalent(which(tr@pvalues < 0.05), 1:2)
  expect_equivalent(which(tr@gof.decimal), 1:4)

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

# bife (bife) ----
test_that("extract bife objects from the bife package", {
  testthat::skip_on_cran()
  skip_if_not_installed("bife", minimum_version = "0.7")
  require("bife")

  set.seed(12345)
  mod <- bife(LFP ~ I(AGE^2) + log(INCH) + KID1 + KID2 + KID3 + factor(TIME) | ID, psid)

  tr <- extract(mod)

  expect_equivalent(sum(tr@coef), 4.61698, tolerance = 1e-2)
  expect_equivalent(sum(tr@se), 2.443803, tolerance = 1e-2)
  expect_equivalent(sum(tr@pvalues), 1.359822, tolerance = 1e-2)
  expect_equivalent(tr@gof, c(-3026.476, 6052.951, 5976), tolerance = 1e-2)
  expect_length(tr@gof.names, 3)
  expect_length(tr@coef, 13)
  expect_equivalent(which(tr@pvalues < 0.05), c(1:4, 8:13))
  expect_equivalent(which(tr@gof.decimal), 1:2)
  expect_equivalent(dim(matrixreg(mod)), c(30, 2))
})

# brmsfit (brms) ----
test_that("extract brmsfit objects from the brms package", {
  testthat::skip_on_cran()
  skip_if_not_installed("brms", minimum_version = "2.8.8")
  skip_if_not_installed("coda", minimum_version = "0.19.2")
  require("brms")
  require("coda")

  # example 2 from brm help page; see ?brm
  sink("/dev/null")
  suppressMessages(fit2 <- brm(rating ~ period + carry + cs(treat),
                               data = inhaler, family = sratio("logit"),
                               prior = set_prior("normal(0,5)"), chains = 1))
  sink()

  suppressWarnings(tr <- extract(fit2))
  expect_length(tr@gof.names, 4)
  expect_length(tr@coef, 8)
  expect_length(tr@se, 8)
  expect_length(tr@pvalues, 0)
  expect_length(tr@ci.low, 8)
  expect_length(tr@ci.up, 8)
  expect_equivalent(which(tr@gof.decimal), c(1, 3, 4))

  # example 1 from brm help page; see ?brm
  bprior1 <- prior(student_t(5,0,10), class = b) + prior(cauchy(0,2), class = sd)
  sink("/dev/null")
  fit1 <- suppressMessages(brm(count ~ zAge + zBase * Trt + (1|patient),
                               data = epilepsy,
                               family = poisson(),
                               prior = bprior1))
  sink()

  suppressMessages(tr <- extract(fit1, use.HDI = TRUE, reloo = TRUE))
  expect_length(tr@gof.names, 5)
  expect_length(tr@coef, 5)
  expect_length(tr@se, 5)
  expect_length(tr@pvalues, 0)
  expect_length(tr@ci.low, 5)
  expect_length(tr@ci.up, 5)
  expect_equivalent(which(tr@gof.decimal), c(1:2, 4:5))
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
  expect_equivalent(tr@coef, c(0.18, 0.43, 0.51), tolerance = 1e-2)
  expect_equivalent(tr@se, c(0.158, 0.053, 0.056), tolerance = 1e-2)
  expect_equivalent(tr@pvalues, c(0.25, 0.00, 0.00), tolerance = 1e-2)
  expect_equivalent(tr@gof, c(0.68, 0.68, 180, 0.04), tolerance = 1e-2)
  expect_length(tr@gof.names, 4)
  expect_length(tr@coef, 3)
  expect_equivalent(which(tr@pvalues < 0.05), 2:3)
  expect_equivalent(which(tr@gof.decimal), c(1, 2, 4))
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

  expect_equivalent(tr@coef, c(1.091, -1.106, 1.123), tolerance = 1e-2)
  expect_equivalent(tr@se, c(0.024, 0.024, 0.024), tolerance = 1e-2)
  expect_equivalent(tr@pvalues, c(0, 0, 0), tolerance = 1e-2)
  expect_equivalent(tr@gof, c(15970.75, 19940.00, 997.00, 20.00), tolerance = 1e-2)
  expect_length(tr@gof.names, 4)
  expect_length(tr@coef, 3)
  expect_equivalent(which(tr@pvalues < 0.05), 1:3)
  expect_equivalent(which(tr@gof.decimal), 1)
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

  # without fixed effects
  OLS1 <- felm(Sepal.Length ~ Sepal.Width |0|0|0, data = iris)
  tr1 <- extract(OLS1)
  expect_length(tr1@gof, 5)
})

# gamlssZadj (gamlss.inf) ----
test_that("extract gamlssZadj objects from the gamlss.inf package", {
  testthat::skip_on_cran()
  skip_if_not_installed("gamlss.inf", minimum_version = "1.0.1")
  require("gamlss.inf")

  set.seed(12345)
  sink("/dev/null")
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
  skip_if_not_installed("lme4")
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
  skip_if_not_installed("lme4")
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

# multinom (nnet), mlogit (mlogit), mnlogit (mnlogit) ----
test_that("extract multinom objects from the nnet package and mlogit and mnlogit objects", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("carData", minimum_version = "3.0.2")
  testthat::skip_if_not_installed("nnet", minimum_version = "7.3.12")
  require("nnet")
  set.seed(12345)
  data(WVS, package = "carData")
  WVS$poverty <- as.factor(as.character(WVS$poverty))
  sink("/dev/null")
  model1 <- multinom(poverty ~ religion + degree + country + age + gender, data = WVS)
  sink()
  tr1 <- extract(model1)

  expect_equivalent(sum(abs(tr1@coef)), 8.387416, tolerance = 1e-2)
  expect_equivalent(sum(tr1@se), 1.624519, tolerance = 1e-2)
  expect_equivalent(sum(tr1@pvalues), 2.961305, tolerance = 1e-2)
  expect_equivalent(sum(tr1@gof), 30608.87, tolerance = 1e-2)
  expect_length(tr1@coef, 16)
  expect_length(tr1@gof, 6)
  expect_equivalent(which(tr1@gof.decimal), 1:4)

  tr1b <- extract(model1, beside = TRUE)
  expect_type(tr1b, "list")
  expect_length(tr1b, 2)
  expect_equivalent(length(tr1b[[1]]@coef), length(tr1b[[1]]@coef))
  expect_length(tr1b[[1]]@coef, 8)

  testthat::skip_if_not_installed("MASS")
  sink("/dev/null")
  example(birthwt, package = "MASS")
  bwt.mu <- multinom(low ~ ., bwt)
  sink()
  tr1c <- extract(bwt.mu)

  expect_equivalent(sum(abs(tr1c@coef)), 8.117209, tolerance = 1e-2)
  expect_equivalent(sum(tr1c@se), 5.314884, tolerance = 1e-2)
  expect_equivalent(sum(tr1c@pvalues), 2.295426, tolerance = 1e-2)
  expect_equivalent(sum(tr1c@gof), 759.248, tolerance = 1e-2)
  expect_length(tr1c@coef, 11)
  expect_length(tr1c@gof, 6)
  expect_equivalent(which(tr1c@gof.decimal), 1:4)
  expect_equivalent(tr1c@gof[6], 2)
  expect_equal(dim(matrixreg(tr1c)), c(29, 2))

  testthat::skip_if_not_installed("mlogit", minimum_version = "1.0.2")
  require("mlogit")
  set.seed(12345)
  WVS2 <- mlogit.data(WVS, choice = "poverty", shape = "wide")
  model2 <- mlogit(poverty ~ 0 | religion + degree + country +
                     age + gender, data = WVS2)
  tr2 <- extract(model2)

  expect_equivalent(sum(abs(tr1@coef)), sum(abs(tr2@coef)), tolerance = 1e-2)
  expect_equivalent(sum(tr1@se), sum(tr2@se), tolerance = 1e-2)
  expect_equivalent(sum(tr1@pvalues), sum(tr2@pvalues), tolerance = 1e-2)
  expect_length(tr2@gof, 4)
  expect_equivalent(which(tr2@gof.decimal), 1:2)

  tr2b <- extract(model2, beside = TRUE)
  expect_type(tr2b, "list")
  expect_length(tr2b, 2)
  expect_equivalent(length(tr2b[[1]]@coef), length(tr2b[[1]]@coef))
  expect_length(tr2b[[1]]@coef, 8)

  tr2c <- extract(model2, include.iterations = TRUE)
  expect_length(tr2c@gof, 6)
  expect_equivalent(which(tr2c@gof.decimal), c(1:2, 6))

  testthat::skip_if_not_installed("mnlogit", minimum_version = "1.2.6")
  require("mnlogit")
  set.seed(12345)
  model3 <- mnlogit(poverty ~ 1 | religion + degree +
                      country + age + gender, data = WVS2)
  tr3 <- extract(model3)

  expect_equivalent(sum(abs(tr3@coef)), sum(abs(tr2@coef)), tolerance = 1e-2)
  expect_equivalent(sum(tr3@se), sum(tr2@se), tolerance = 1e-2)
  expect_equivalent(sum(tr3@pvalues), sum(tr2@pvalues), tolerance = 1e-2)
  expect_length(tr3@gof, 4)
  expect_equivalent(which(tr3@gof.decimal), 1:2)

  tr3b <- extract(model3, beside = TRUE)
  expect_type(tr3b, "list")
  expect_length(tr3b, 2)
  expect_equivalent(length(tr3b[[1]]@coef), length(tr3b[[1]]@coef))
  expect_length(tr3b[[1]]@coef, 8)

  tr3c <- extract(model3, include.iterations = TRUE)
  expect_length(tr3c@gof, 7)
  expect_equivalent(which(tr3c@gof.decimal), c(1:2, 6:7))

  mr1 <- matrixreg(list(model1, model2, model3))
  expect_equivalent(dim(mr1), c(103, 4))
  mr2 <- matrixreg(list(model1, model2, model3), beside = TRUE)
  expect_equivalent(dim(mr2), c(25, 7))
})

# nlmerMod (lme4) ----
test_that("extract nlmerMod objects from the lme4 package", {
  testthat::skip_on_cran()
  skip_if_not_installed("lme4")
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
