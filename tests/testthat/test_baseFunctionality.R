context("Test base functionality")

n <- 100
set.seed(123)
dt <- data.frame(y = rnorm(n),
                 fe1 = sample(x = letters[1:10], size = n, replace = T),
                 fe2 = sample(x = letters[11:12], size = n, replace = T),
                 x1 = runif(n),
                 x2 = runif(n) + 4)

res.x1 <- lm(y ~ x1, data = dt)
res.x2 <- lm(y ~ x2, data = dt)
res.x1x2 <- lm(y ~ x1 + x2, data = dt)
res.noInt <- lm(y ~ - 1 + x1 + x2, data = dt)

texreg(list(res.x1, res.x2), file = "test_tab1.tex")
texreg(list(res.x1, res.x2), file = "test_tab2.tex", booktabs = TRUE)
texreg(list(res.x1, res.x2), file = "test_tab3.tex", booktabs = TRUE,
       dcolumn = TRUE)
texreg(list(res.x1, res.x2), file = "test_tab4.tex", single.row = TRUE,
       booktabs = TRUE, dcolumn = TRUE)
texreg(list(res.x1, res.x2), file = "test_tab5.tex", 
       booktabs = TRUE, dcolumn = TRUE,
       custom.model.names = c("Test1", "Test2"))


##
## test all test tables from above against the expected tables
test_that("Test equality of output files", {
    for(i in 1:5) {
        expect_equal(readLines(paste0("test_tab", i, ".tex")),
                     readLines(paste0("expect_tab", i, ".tex")))
    }
})
