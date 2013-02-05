## some calculations that LLA has to replicate, done in R
##
## this code is not run for the Lisp tests

## least squares
x <- matrix(c(23,23,22,21,25,20,29,32,24,29),5,2,byrow=TRUE)
y <- x %*% c(1,2) + c(-2,-1,0,1,2)
fit <- lm.fit(x,y)
sum(fit$residuals^2)
fit.full <- lm(y~x-1)


## cholesky
a <- matrix(c(1,2,0,3),2,2,byrow=TRUE)
b <- t(a) %*% a
chol(t(a) %*% a)


## eigenvalues -- hermitian
a <- matrix(1:4,2,2,byrow=TRUE)
aa <- t(a) %*% a
eigen(aa)

## svd
La.svd(matrix(1:6,3,2,byrow=TRUE))
