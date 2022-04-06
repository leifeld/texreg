context("praise function")
suppressPackageStartupMessages(library("texreg"))

test_that("praise gives useful error messages if the fields were not filled out properly", {
  expect_error(praise(organization = "test"), "argument \"academic_user\" is missing, with no default")
  expect_error(praise(academic_user = TRUE), "argument \"organization\" is missing, with no default")
  expect_error(praise(academic_user = TRUE, organization = FALSE), "'organization' must be the name of your organization.")
  expect_error(praise(academic_user = 24), "'academic_user' must be TRUE if you are at a university or research institute and FALSE otherwise.")
  expect_error(praise(academic_user = TRUE, organization = "test", name = 25, general_praise = "test"), "Optional 'name' argument:")
  expect_error(praise(academic_user = TRUE, organization = "test"),
               "At least one of the arguments 'general_praise', 'increase_productivity', or 'increase_quality' must be provided")
  expect_error(praise(academic_user = TRUE, organization = "test", general_praise = NA),
               "'general_praise' is a free-form text argument where you can provide any kind of praise.")
  expect_error(praise(academic_user = TRUE, organization = "test", increase_productivity = NA),
               "'increase_productivity' argument: Tell us how texreg makes you more productive.")
  expect_error(praise(academic_user = TRUE, organization = "test", increase_quality = NA),
               "'increase_quality' argument: Tell us how texreg increases the quality of your work or reporting.")
  expect_error(praise(academic_user = TRUE, organization = "test", general_praise = "test", start_using = TRUE),
               "Optional 'start_using' argument: Tell us when you started using texreg")
  expect_error(praise(academic_user = TRUE, organization = "test", general_praise = "test", where_learn = 25),
               "Optional 'where_learn' argument: Tell us where or how you learned about texreg")
  expect_error(praise(academic_user = TRUE, organization = "test", general_praise = "test", contact_details = 25),
               "Optional 'contact_details' argument: How can we reach you if we need to ask for further testimony for our reporting purposes?")
  expect_error(praise(academic_user = TRUE, organization = "test", general_praise = "test", models = 25),
               "Optional 'models' argument: What kinds of statistical models do you use in your work most often?")
  expect_error(praise(academic_user = TRUE, organization = "test", general_praise = "test", num_users = TRUE),
               "Optional 'num_users' argument: How many other texreg users do you know approximately[?] Supply a number or some text.")
})