context("Test extract.alpaca")

countries <- letters[1:10]
years <- 2000:2010

dta <- expand.grid(exp = countries,
                   imp = countries,
                   year = years)
dta$tradeflow <- rpois(nrow(dta), lambda = sqrt(nrow(dta)))
dta$tariffs <- abs(rnorm(nrow(dta)))
dta$love <- abs(rnorm(nrow(dta)))
dta$expFE <- with(dta, interaction(exp, year))
dta$impFE <- with(dta, interaction(imp, year))
dta$yearFE <- factor(dta$year)
dta$bilFE <- with(dta, interaction(exp, imp))

model <- feglm(tradeflow ~ tariffs + love | expFE + impFE + yearFE + bilFE,
             data = dta, family = poisson())

s <- summary(model)
