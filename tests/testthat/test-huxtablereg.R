context("huxtablereg")
suppressPackageStartupMessages(library("texreg"))

skip_if_not_installed("huxtable")

library("huxtable")

m1 <- lm(mpg ~ gear, mtcars)
m2 <- lm(mpg ~ gear + cyl, mtcars)

test_that("huxtablereg gives useful error message if huxtable not installed", {
  skip_if_not_installed("huxtable", minimum_version = "5.0.0")
  with_mock(requireNamespace = function (...) return(FALSE), {
    expect_error(huxtablereg(list(m1, m2)), regex = "huxtable")
  })
})

test_that("huxtablereg works", {
  skip_if_not_installed("huxtable", minimum_version = "5.0.0")
  expect_equivalent(dim(as.data.frame(huxtablereg(list(m1, m2)))), c(10, 3))
})

test_that("huxtablereg works with beside = TRUE", {
  skip_if_not_installed("systemfit")

  library(systemfit)
  data("Kmenta")
  eqDemand <- consump ~ price + income
  eqSupply <- consump ~ price + farmPrice + trend
  system <- list(demand = eqDemand, supply = eqSupply)
  fitols <- systemfit(system, data = Kmenta)

  h <- huxtablereg(fitols, beside = TRUE)
  expect_equal(ncol(h), 3)
})
