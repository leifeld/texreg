context("plotreg")
suppressPackageStartupMessages(library("texreg"))

skip_if_not_installed("ggplot2")
require("ggplot2")

test_that("plotreg works", {
  skip_if_not_installed("ggplot2", minimum_version = "3.3.1")
  # linear model p-values
  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight <- c(ctl, trt)
  m1 <- lm(weight ~ group)
  m2 <- lm(weight ~ group - 1)
  p1 <- plotreg(list(m1, m2))
  p2 <- plotreg(list(m1, m2), ci.force = TRUE)
  p3 <- plotreg(list(m1, m2), custom.title = "My plot")
  p4 <- plotreg(list(m1, m2), custom.note = "My note")

  expect_true(ggplot2::is.ggplot(p1))
  expect_true(ggplot2::is.ggplot(p2))
  expect_true(ggplot2::is.ggplot(p3))
  expect_true(ggplot2::is.ggplot(p4))

  expect_s3_class(p1$data, "data.frame")
  expect_s3_class(p1, "gg")
  expect_equal(dim(p1$data), dim(readRDS("../files/PlotDataFrame.RDS"))) # test if data frame is correctly constructed
   #saveRDS(p1$data, "../files/PlotDataFrame.RDS")

  expect_identical(p1$labels$y, "Bars denote SEs (95%). Circle points denote significance.")# tests if it is printing the correct xlab
  expect_identical(p2$labels$y, "Bars denote CIs (95%) computed from SEs. Circle points denote significance.") # tests if it is printing the correct xlab
  expect_identical(p3$labels$title, "My plot") # test if title is printed
  expect_identical(p4$labels$y, "My note") # test if note is printed
})

test_that("plotreg works with confidence intervals using the biglm package", {
  skip_if_not_installed("ggplot2", minimum_version = "3.3.1")
  skip_if_not_installed("biglm")
  require("biglm")
  skip_if_not_installed("DBI")
  data("trees")
  ff <- log(Volume) ~ log(Girth) + log(Height)
  a <- biglm::bigglm(ff, data = trees, chunksize = 10, sandwich = TRUE)
  gg <- log(Volume) ~ log(Girth) + log(Height) + offset(2 * log(Girth) + log(Height))
  b <- biglm::bigglm(gg, data = trees, chunksize = 10, sandwich = TRUE)
  p6 <- plotreg(list(a, b))
  expect_true(ggplot2::is.ggplot(p6))
  expect_equal(dim(p6$data), dim(readRDS("../files/PlotDataFrameCI.RDS")))
  # saveRDS(p6$data, "../files/PlotDataFrameCI.RDS")
  expect_identical(p6$labels$y, "Bars denote SEs (95%). Circle points denote significance.")
})

test_that("plotreg -odds ratio", {
  skip_if_not_installed("ggplot2", minimum_version = "3.3.1")
  # linear model p-values
  # the model is not actually suitable for odds ratio estimation, but it
  # makes no difference for testing plotreg
  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight <- c(ctl, trt)
  m1 <- lm(weight ~ group)
  p7 <- plotreg(m1, override.coef = exp(m1$coefficients), ci.force = TRUE, ci.test = 1)

  expect_true(ggplot2::is.ggplot(p7))
  expect_s3_class(p7$data, "data.frame")
  expect_s3_class(p7, "gg")
  expect_equal(dim(p7$data), dim(readRDS("../files/PlotDataFrameOR.RDS"))) # test if data frame is correctly constructed
  #saveRDS(p1$data, "../files/PlotDataFrameOR.RDS")
  class(p7$data$signif)
  # tests if significance is computed according to the value provided with ci.test
  expect_identical(as.character(p7$data$signif), as.character(c(TRUE, FALSE)))
})