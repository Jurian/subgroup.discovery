
library(subgroup.discovery)


data(pima)
p.cov <- prim.cover(class ~ ., data = pima, peeling.quantile = 0.05, min.support = 0.05, plot = TRUE)
summary(p.cov)
plot(p.cov)



data(pima)
p.div <- prim.diversify(class ~ ., data = pima, n = 10, peeling.quantile = 0.05, min.support = 0.05, plot = TRUE)
summary(p.div)
plot(p.div)




data(ames)
p.cov <- prim.cover(SalePrice ~ . - PID - Order, ames, peeling.quantile = 0.05, min.support = 0.1, plot = TRUE)
summary(p.cov)
plot(p.cov)




data(ames)
p.div <- prim.diversify(SalePrice ~ . - PID - Order, ames, n = 10, peeling.quantile = 0.05, min.support = 0.05, plot = TRUE)
summary(p.div)
plot(p.div)

