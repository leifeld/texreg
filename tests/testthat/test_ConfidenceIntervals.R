context("Test confidence intervals")

data(cars)
model <- lm(dist ~ speed, data=cars)
out <- texreg(model, ci.force=T)

## very helpful:
## sprintf(out)

expect_true(grepl("\\(Intercept\\) *& *\\$-17\\.58\\^\\{\\*\\}\\$ *\\\\ *\n",
                  out, perl = TRUE))
expect_true(grepl("\\(Intercept\\) *& *\\$-17\\.58\\^\\{\\*\\}\\$ *\\\\ *& *\\$[-30.83;\\ -4.33]\\$", out))
expect_true(grepl("\\(Intercept\\) *& *\\$-17\\.58\\^\\{\\*\\}\\$ *\\\\ *& *\\$[-30.83;\\ -4.33]\\$ *\\\\", out))

