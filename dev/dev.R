
library(subgroup.discovery)

data(pima)
pima.sample <- sample(nrow(pima), 0.75*nrow(pima))
pima <- prim.data.prepare(pima)
pima.peel <- prim(class ~ ., data = pima[pima.sample,], peeling.quantile = 0.3, min.support = 0.4)
summary(pima.peel)
plot(pima.peel)
pima.predict <- predict(pima.peel, pima[-pima.sample,])
summary(pima.predict)
plot(pima.predict)

data(ames)
ames.sample <- sample(nrow(ames), 0.75*nrow(ames))
ames <- prim.data.prepare(ames)
ames.peel <- prim(SalePrice ~ . - PID - Order, data = ames[ames.sample,], peeling.quantile = 0.05, min.support = 0.1)
summary(ames.peel)
plot(ames.peel)
ames.predict <- predict(ames.peel, ames[-ames.sample,])
summary(ames.predict)
plot(ames.predict)

