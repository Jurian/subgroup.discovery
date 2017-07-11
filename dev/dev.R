
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
p.cov <- prim.cover(X = ames[,-c(1,2,82)], y = ames$SalePrice, peeling.quantile = 0.05, min.support = 0.1, plot = TRUE)
summary(p.cov)
plot(p.cov)




data(ames)
p.div <- prim.diversify(X = ames[,-c(1,2,82)], y = ames$SalePrice, n = 10, peeling.quantile = 0.05, min.support = 0.05, plot = TRUE)
summary(p.div)
plot(p.div)



#X <- stats::model.frame(formula = SalePrice ~ . - PID - Order, data = ames)
