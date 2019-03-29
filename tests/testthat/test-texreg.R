context("Test of various arguments in texreg")
library("texreg")

data("iris")
model1 <- lm(Sepal.Width ~ Petal.Width, data = iris)

test_that("no stars, single.row, and dcolumn work together", {
  expect_match(texreg(model1,
                      single.row = TRUE,
                      dcolumn = TRUE,
                      stars = numeric()),
               regexp = "l D{\\)}{\\)}{11\\)0}}",
               perl = TRUE)
})
