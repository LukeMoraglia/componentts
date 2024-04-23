set.seed(42)
X <- matrix(sample.int(20, 1000*100, replace = TRUE), 1000, 100)
set.seed(69)
Y <- matrix(sample.int(20, 1000*50, replace = TRUE), 1000, 50)
X_s <- scale_SS1(X)
Y_s <- scale_SS1(Y)

res1 <- t(X_s) %*% Y_s
res2 <- componentts:::mat_mult(t(X_s), Y_s)

test_that("fast mat mult", {
  expect_equal(res2, res1)
})

tX_s <- t(X_s)
bench::mark(
  tX_s %*% Y_s,
  componentts:::mat_mult(tX_s, Y_s),
  Rfast::mat.mult(tX_s, Y_s)
)
# 7/19/23: 3.9ms, 1ms, 7.8ms
