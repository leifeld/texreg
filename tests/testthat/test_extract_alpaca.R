context("Test extract.alpaca")

library(alpaca)

set.seed(1234)
countries <- letters[1:10]
years <- 2000:2010

dta <- expand.grid(exp = countries,
                   imp = countries,
                   year = years)
dta$tradeflow <- rpois(nrow(dta), lambda = sqrt(nrow(dta)))
dta$tariffs <- abs(rnorm(nrow(dta)))
dta$var1 <- abs(rnorm(nrow(dta)))
dta$var2 <- abs(rnorm(nrow(dta)))
dta$expFE <- with(dta, interaction(exp, year))
dta$impFE <- with(dta, interaction(imp, year))
dta$yearFE <- factor(dta$year)
dta$bilFE <- with(dta, interaction(exp, imp))

mod1 <- feglm(tradeflow ~ tariffs + var1 | expFE + impFE + yearFE + bilFE,
              data = dta, family = poisson())
mod2 <- feglm(tradeflow ~ tariffs + var2 | expFE + impFE + yearFE + bilFE,
              data = dta, family = poisson())


test_that("Test extract.alpaca function", {
    res <- texreg::extract.alpaca(mod1)

    ## Number of slots
    expect_equal(length(slotNames(res)), 10)

    ## Names of coefficients
    expect_equal(res@coef.names, c("tariffs", "var1"))

    ## Values of coefficients
    expect_equal(res@coef, c(tariffs = -0.00690423707039415,
                             var1 = 0.00179111096513545),
                 tol = 1e-8)

    ## Values of p-values
    expect_equal(res@pvalues, c(tariffs = 0.487814242755162,
                                var1 = 0.861041947556746),
                 tol = 1e-8)

    ## Value of deviance, number of observations
    expect_equal(res@gof, c(815.123394, 1100), tol = 1e-8)
})

test_that("Table has correct dimensions, three cols, five rows", {
    res <- texreg(list(mod1, mod2))

    ## three columns
    expect_true(grepl("\\\\begin\\{tabular\\}\\{l c c \\}", res))
    expect_true(grepl("\\\\multicolumn\\{3\\}\\{l\\}", res))

    ## five rows
    expect_true(grepl("\\\ntariffs", res))
    expect_true(grepl("\\\nvar1", res))
    expect_true(grepl("\\\nvar2", res))
    expect_true(grepl("\\\nDeviance", res))
    expect_true(grepl("\\\nNum. obs.", res))
})
