context("plotreg")
suppressPackageStartupMessages(library("texreg"))

skip_if_not_installed("ggplot2")
require("ggplot2")

test_that("plotreg works", {
  # linear model p-values
  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight <- c(ctl, trt)
  m1 <- lm(weight ~ group)
  m2 <- lm(weight ~ group - 1)
  p1 <- plotreg(list(m1,m2))
  p2 <- plotreg(list(m1,m2), use.se = TRUE)
  p3 <- plotreg(list(m1,m2), custom.title = "My Plot")
  p4 <- plotreg(list(m1,m2), custom.note = "My note")

  expect_true(ggplot2::is.ggplot(p1))
  expect_true(ggplot2::is.ggplot(p2))
  expect_true(ggplot2::is.ggplot(p3))
  expect_true(ggplot2::is.ggplot(p4))

  expect_equal(PlotDataFrame <- p1$data, readRDS("../files/PlotDataFrame.RDS")) # test if data frame is correctly constructed
  # saveRDS(PlotDataFrame, "../files/PlotDataFrame.RDS")

  expect_identical(p1$labels$y, "Bars denote CIs. Circle points denote significance.") # tests if it is printing the correct xlab
  expect_identical(p2$labels$y, "Bars denote SEs. Circle points denote significance.") # tests if it is printing the correct xlab
  expect_identical(p3$labels$title, "My Plot") # test if title is printed
  expect_identical(p4$labels$y, "My note") # test if note is printed

  p5 <- plotreg(list(m1, m2), ci.inner = 0)
  expect_equal(p5$data$lower.inner, p5$data$upper.inner, tolerance = 1e-3) # test if inner CIs can be switched off
})

test_that("plotreg works with confidence intervals using the biglm package", {
  skip_if_not_installed("biglm")
  require("biglm")
  skip_if_not_installed("DBI")
  data("trees")
  ff <- log(Volume) ~ log(Girth) + log(Height)
  a <- biglm::bigglm(ff, data = trees, chunksize = 10, sandwich = TRUE)
  gg <- log(Volume) ~ log(Girth) + log(Height) + offset(2 * log(Girth) + log(Height))
  b <- biglm::bigglm(gg, data = trees, chunksize = 10, sandwich = TRUE)
  p6 <- plotreg(list(a, b))
  expect_true(is.ggplot(p6))
  expect_equal(PlotDataFrameCI <- p6$data, readRDS("../files/PlotDataFrameCI.RDS"))
  # saveRDS(PlotDataFrameCI, "../files/PlotDataFrameCI.RDS")
  expect_identical(p6$labels$y, "Bars denote CIs. Circle points denote significance.")
})