#library(bench)

set.seed(42)
X <- matrix(sample.int(20, 20*3, replace = TRUE), 20, 3)
set.seed(69)
Y <- matrix(sample.int(20, 20*4, replace = TRUE), 20, 4)
X_s <- scale_SS1(X)
Y_s <- scale_SS1(Y)

Rx <- t(X_s) %*% X_s

bench::mark(
  GSVD::invsqrt_psd_matrix(Rx),
  fast_invsqrt_psd_matrix(Rx, symmetric = TRUE)
)

# X2 <- cbind(X_s, X_s[,1])
# Rx <- t(X2) %*% X2
# res <- GSVD::tolerance_eigen(Rx)
# invsqrt <- t(t(res$vectors) * (1/sqrt(res$values)) ) %*% t(res$vectors)
# inv <- invsqrt %*% invsqrt
# inv
# Rx %*% inv
#
#
# res <- eigen(Rx)
# invsqrt <- t(t(res$vectors) * (1/sqrt(res$values)) ) %*% t(res$vectors)
# inv <- invsqrt %*% invsqrt
# inv
# Rx %*% inv

# Using eigen vs tolerance_eigen gives different inverses, depending on if the matrix is singular or not
