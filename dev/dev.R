
library(subgroup.discovery)


data(pima)
p.cov <- prim.cover(class ~ ., data = pima, plot = TRUE, optimal.box = "2se")
summary(p.cov)
plot(p.cov)



data(pima)
p.div <- prim.diversify(class ~ ., data = pima, n = 100, plot = TRUE, parallel = TRUE, optimal.box = "2se")
summary(p.div)
plot(p.div)




data(ames)
p.cov <- prim.cover(SalePrice ~ . - PID - Order, ames, plot = TRUE)
summary(p.cov)
plot(p.cov)




data(ames)
p.div <- prim.diversify(SalePrice ~ . - PID - Order, ames, n = 3, plot = TRUE)
summary(p.div)
plot(p.div)

