# dataframe

ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
weight <- c(ctl, trt)
m1 <- lm(weight ~ group)
m2 <- lm(weight ~ group - 1)
p <- plotreg(list(m1,m2))
#saveRDS(PlotDataFrame, "../files/PlotDataFrame.RDS)
#myplot <- plotreg(list(m1,m2))
#PlotDataFrame <- myplot$data
#saveRDS(PlotDataFrame, "/Users/claudiazucca/Documents/texreg/tests/files/PlotDataFrame.RDS")

test_that("Dataframe is correctly constructed",{
    #p <- plotreg(list(m1,m2))
    expect_true(is.ggplot(p))
    expect_equal(p$data, readRDS("../files/PlotDataFrame.RDS"))
})

p$data

#######################################################


# library(ggplot2)
# library(scales) # for percent()
# library(testthat)
# 
# # First, 'correct' data frame
# df <- data.frame(
#     Response   = LETTERS[1:5],
#     Proportion = c(0.1,0.2,0.1,0.2,0.4)
# )
# 
# # Second data frame where column has 'wrong' name that does not match aes()
# df2 <- data.frame(
#     x          = LETTERS[1:5],
#     Proportion = c(0.1,0.2,0.1,0.2,0.4)
# )
# 
# plot_fun <- function(df) {
#     p1 <- ggplot(df, aes(Response, Proportion)) +
#         geom_bar(stat='identity') + 
#         scale_y_continuous(labels = percent)
#     return(p1)
# }
# 
# # All tests succeed
# test_that("Scale is labelled 'Proportion'",{
#     p <- plot_fun(df)
#     expect_true(is.ggplot(p))
#     expect_identical(p$labels$y, "Proportion")
#     
#     p <- plot_fun(df2)
#     expect_true(is.ggplot(p))
#     expect_identical(p$labels$y, "Proportion")
# })
# 
# # Second test with data frame df2 fails
# test_that("Printing ggplot object actually works",{
#     p <- plot_fun(df)
#     expect_error(print(p), NA)
#     
#     p <- plot_fun(df2)
#     expect_error(print(p), NA)
# })
# #> Error: Test failed: 'Printing ggplot object actually works'
# #> * `print(p)` threw an error.
# #> Message: object 'Response' not found
# #> Class:   simpleError/error/condition