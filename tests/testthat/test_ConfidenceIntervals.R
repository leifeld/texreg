context("Test confidence intervals")

data(cars)
model <- lm(dist ~ speed, data=cars)
out <- texreg(model, ci.force=T)

## The important difference of the next two tests is the ending:
## In the first test, there is the beginning of an additional CI
## block [-30.83; ...] (which, of course, should not be there).
## In the second test, there is the -- correct -- double
## backslash and the newline character.

## Is currently (as of 1.36.24) true, but should be false
expect_false(grepl("\\(Intercept\\) *& *\\$-17\\.58\\^\\{\\*\\}\\$ *\\\\\\\\\n *& *\\$\\[-30.83;\\\\ -4.33\\]\\$ *\\[\\$\\[-30\\.83", out))

## Is currently (as of 1.36.24) false, but should be true
expect_true(grepl("\\(Intercept\\) *& *\\$-17\\.58\\^\\{\\*\\}\\$ *\\\\\\\\\n *& *\\$\\[-30.83;\\\\ -4.33\\]\\$ *\\\\\\\\\n", out))

