# testthat for plotreg

# linear model -pvalues
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
weight <- c(ctl, trt)
m1 <- lm(weight ~ group)
m2 <- lm(weight ~ group - 1)
p <- plotreg(list(m1,m2))
p2 <- plotreg(list(m1,m2), use.se = TRUE)
# tests if dataframe is correctly made
test_that("Dataframe is correctly constructed",{
    expect_true(is.ggplot(p))
    expect_equal(p$data, readRDS("../files/PlotDataFrame.RDS"))
})

# tests if it is printing the correct xlab
test_that("xlab is labelled 'Bars denote CIs. Circle points denote significance.'",{
    expect_true(is.ggplot(p))
    expect_identical(p$labels$y, "Bars denote CIs. Circle points denote significance.")
})

# tests if it is printing the correct xlab
test_that("xlab is labelled 'Bars denote SEs. Circle points denote significance.'",{
    expect_true(is.ggplot(p2))
    expect_identical(p2$labels$y, "Bars denote SEs. Circle points denote significance.")
})

# biglm model -confidence interval 
skip_if_not_installed("biglm")
require("biglm")
data(trees)
ff <- log(Volume) ~ log(Girth) + log(Height)
a <- bigglm(ff, data = trees, chunksize = 10, sandwich = TRUE)
gg <- log(Volume) ~ log(Girth) + log(Height) + offset(2 * log(Girth) + log(Height))
b <- bigglm(gg, data = trees, chunksize = 10, sandwich = TRUE)
p3 <- plotreg(list(a, b))

# tests if dataframe is correctly made
test_that("biglm tests",{
    expect_true(is.ggplot(p3))
    expect_equal(p3$data, readRDS("../files/PlotDataFrameCI.RDS"))
    expect_identical(p3$labels$y, "Bars denote CIs. Circle points denote significance.")
})

# # tests if it is printing the correct xlab
# test_that("xlab is labelled 'Bars denote CIs. Circle points denote significance.'",{
#     expect_true(is.ggplot(p3))
#     expect_identical(p3$labels$y, "Bars denote CIs. Circle points denote significance.")
# })








