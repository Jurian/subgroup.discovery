
library(bump.hunting)
library(data.table)

dat <- fread("inst/extdata/credit.csv")
dat$gender <- factor(dat$gender)

X <- dat[,-6]
y <- dat$class


p.train <- prim.default(X = X, y = y, peeling.quantile  = 0.01, min.support = 0.4)
plot(p.train)











dat <- fread("inst/extdata/pima.csv")
X <- dat[,-9]
y <- dat$class


train <- sample(1:nrow(X), nrow(X) * 0.6)

#p.train <- prim.formula(class ~ ., dat[train], 0.01, 0.4)
p.train <- prim.default(X = X[train], y = y[train], peeling.quantile = 0.01, min.support = 0.4)
p.test <- prim.test(p.train, X[-train], y[-train])

plot(p.train)
plot(p.test)


t <- prim.cover(X, y, 0.01, 0.4, 0.6)

