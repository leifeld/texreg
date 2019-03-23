context("huxtablereg")

skip_if_not_installed("huxtable")

library("texreg")
library("huxtable")

m1 <- lm(mpg ~ gear, mtcars)
m2 <- lm(mpg ~ gear + cyl, mtcars)

test_that("huxtablereg gives useful error message if huxtable not installed", {
  with_mock(requireNamespace = function (...) return(FALSE), {
    expect_error(huxtablereg(list(m1, m2)), regex = "huxtable")
  })
})

# # WORKS
# test_that("screenreg works", {
#   expect_equal(output <- screenreg(list(m1, m2)),
#                readRDS("../files/screenreg2.RDS"))
# })
# # saveRDS(output, "../files/screenreg2.RDS")
#
# # FAILS
# test_that("huxtablereg works", {
#   expect_equal(output <- huxtablereg(list(m1, m2)),
#                readRDS("../files/huxtablereg.RDS"))
# })
# # saveRDS(output, "../files/huxtablereg.RDS")
