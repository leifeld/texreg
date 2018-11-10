context("Test base functionality")

## generate some test random data
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

## test for correct default coefficient names and values
out <- texreg(res.x1)
expect_true(grepl("\\(Intercept\\) *& *\\$0.17\\$ *\\\\", out))
expect_true(grepl("x1 +& *\\$-0.15\\$", out))
expect_true(grepl("& Model 1", out))

## test for correct customisation
out <- texreg(res.x1, custom.model.names = "Custom 1",
              custom.coef.names = c("someName", "someOther"))
expect_true(grepl("& Custom 1", out))
expect_true(grepl("someName *& *\\$0.17\\$ *\\\\", out))
expect_true(grepl("someOther +& *\\$-0.15\\$", out))

