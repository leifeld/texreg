
# Helpful commands in R for writing tests with regexes

## sprintf

Prints the string as it is, included the newline characters

```r
model <- lm(dist ~ speed, data=cars)
out <- texreg(model, ci.force=T)
sprintf(out)
```

## stringr

`str_extract` prints the strings that are matched with the given
regex. Can be used to verify that you match what you want to match.

```r
library(stringr)
str_extract(out, "*& *\\$-17\\.58\\^\\{\\*\\}\\$ *\\\\\\\\\n *& *\\$\\[-30.83;\\\\ -4.33\\]\\$")

```


