# To get accurate benchmarks with compiled code, do not use devtools::load_all()
# Restart R, reinstall package, restart R, and load package with library(componentts)

set.seed(42)
X <- matrix(sample.int(20, 1000*100, replace = TRUE), 1000, 100)
set.seed(69)
Y <- matrix(sample.int(20, 1000*50, replace = TRUE), 1000, 50)
X_s <- scale_SS1(X)
Y_s <- scale_SS1(Y)

Rx_invsqrt <- fast_invsqrt_psd_matrix(crossprod(X_s), symmetric = TRUE)
Ry_invsqrt <- fast_invsqrt_psd_matrix(crossprod(Y_s), symmetric = TRUE)

matrices <- list(Rx_invsqrt, t(X_s), Y_s, Ry_invsqrt)
split <- componentts:::get_split_mat(matrices)

res1 <- Rx_invsqrt %*% t(X_s) %*% Y_s %*% Ry_invsqrt
res2 <- matrix_chain_multiplication(matrices)
res3 <- matrix_chain_multiplication(matrices, split)

test_that("mat chain mult", {
  expect_equal(res2, res1)
  expect_equal(res3, res1)
})

# CCA situation
bench::mark(
  Rx_invsqrt %*% t(X_s) %*% Y_s %*% Ry_invsqrt,
  matrix_chain_multiplication(matrices),
  matrix_chain_multiplication(matrices, split)
)
# 7/19/23 10.3ms, 3.8ms, 3.6ms
# 7/19/23 12.2ms, 1.4ms, 1.3ms (after adding mat_mult to matrix_chain_multiplication)

# RDA situation
bench::mark(
  Rx_invsqrt %*% t(X_s) %*% Y_s,
  matrix_chain_multiplication(list(Rx_invsqrt, t(X_s), Y_s))
)
# 7/19/23 11.5ms, 5.2ms
# 7/19/23 10.8ms, 1.4ms (after adding mat_mult to matrix_chain_multiplication)

