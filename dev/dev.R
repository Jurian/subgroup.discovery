
library(subgroup.discovery)

data(credit)

y <- credit$class
credit <- credit[,-6]

p.train <- prim.default(X = credit, y = y, peeling.quantile  = 0.1, min.support = 0.4)
plot(p.train)



data(pima)
y <- pima$class
pima <- pima[,-9]

train <- sample(1:nrow(pima), nrow(pima) * 0.75)

#p.train <- prim.formula(class ~ ., dat[train], 0.01, 0.4)
p.train <- prim.default(X = pima[train,], y = y[train], peeling.quantile = 0.01, min.support = 0.01)
p.test <- prim.validate(p.train, pima[-train,], y[-train])

plot(p.train)
plot(p.test)


t <- prim.cover(X, y, 0.01, 0.4, 0.6)



data(pima)

p.cov <- prim.cover(class ~ ., pima, peeling.quantile = 0.05, min.support = 0.05)





data(pima)
p.div <- prim.diversify(X = pima[,-9], y = pima$class, n = 10, peeling.quantile = 0.05, min.support = 0.05)
