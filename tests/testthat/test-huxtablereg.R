context("huxtablereg")

skip_if_not_installed("huxtable")

m1 <- lm(mpg ~ gear, mtcars)
m2 <- lm(mpg ~ gear + cyl, mtcars)

test_that("huxtablereg gives useful error message if huxtable not installed", {
  with_mock(requireNamespace = function (...) return(FALSE), {
    expect_error(huxtablereg(list(m1, m2)), regex = "huxtable")
  })
})


test_that("huxtablereg works", {
  m1 <- lm(mpg ~ gear, mtcars)
  m2 <- lm(mpg ~ gear + cyl, mtcars)
  expect_silent(huxtablereg(list(m1, m2)))
})

