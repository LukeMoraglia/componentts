seed <- sample.int(1000, 1)
print(paste("Random seed:", seed))
set.seed(seed)
row <- 20 #sample(2:100, 1)
col1 <- 10 #sample(2:100, 1)
col2 <- 6 #sample(2:100, 1)
randmat1 <- matrix(rnorm(row * col1), row, col1)
randmat2 <- matrix(runif(row * col2), row, col2)
rmat1_s <- scale_SS1(randmat1)
rmat2_s <- scale_SS1(randmat2)

res_plsc_12 <- plsc(rmat1_s, rmat2_s)
res_plsc_21 <- plsc(rmat2_s, rmat1_s)

test_that("plsc results are equivalent regardless of which matrix is first", {
  expect_equal(res_plsc_12$d, res_plsc_21$d)
  expect_equal(res_plsc_12$l, res_plsc_21$l)
  expect_equal(abs(res_plsc_12$u), abs(res_plsc_21$v))
  expect_equal(abs(res_plsc_12$v), abs(res_plsc_21$u))
  expect_equal(abs(res_plsc_12$lx), abs(res_plsc_21$ly))
  expect_equal(abs(res_plsc_12$ly), abs(res_plsc_21$lx))
  expect_equal(abs(res_plsc_12$p), abs(res_plsc_21$q))
  expect_equal(abs(res_plsc_12$q), abs(res_plsc_21$p))
})

test_that("plsc crossproduct was decomposed", {
  expect_equal(t(rmat1_s) %*% rmat2_s,
               res_plsc_12$p %*% diag(res_plsc_12$d) %*% t(res_plsc_12$q))
})


res_cca <- cca(rmat1_s, rmat2_s)

test_that("cca", {
  expect_equal(t(rmat1_s) %*% rmat2_s,
               res_cca$p %*% diag(res_cca$d) %*% t(res_cca$q))
  expect_equal(t(res_cca$p) %*% res_cca$m %*% res_cca$p,
               diag(1, ncol(res_cca$p), ncol(res_cca$p)))
  expect_equal(t(res_cca$q) %*% res_cca$w %*% res_cca$q,
               diag(1, ncol(res_cca$q), ncol(res_cca$q)))
  expect_equal(diag(t(res_cca$lx) %*% res_cca$ly),
               res_cca$d)
  expect_equal(res_cca$lx,
               rmat1_s %*% res_cca$f)
  expect_equal(res_cca$ly,
               rmat2_s %*% res_cca$g)
  expect_equal(res_cca$f,
               res_cca$m %*% res_cca$p)
  expect_equal(res_cca$g,
               res_cca$w %*% res_cca$q)
})

res_rda <- rda(rmat1_s, rmat2_s)

test_that("rda", {
  expect_equal(t(rmat1_s) %*% rmat2_s,
               res_rda$p %*% diag(res_rda$d) %*% t(res_rda$q))
  expect_equal(t(res_rda$p) %*% res_rda$m %*% res_rda$p,
               diag(1, ncol(res_rda$p), ncol(res_rda$p)))
  expect_equal(t(res_rda$q) %*% res_rda$w %*% res_rda$q,
               diag(1, ncol(res_rda$q), ncol(res_rda$q)))
  expect_equal(diag(t(res_rda$lx) %*% res_rda$ly),
               res_rda$d)
  expect_equal(res_rda$lx,
               rmat1_s %*% res_rda$f)
  expect_equal(res_rda$ly,
               rmat2_s %*% res_rda$g)
  expect_equal(res_rda$f,
               res_rda$m %*% res_rda$p)
  expect_equal(res_rda$g,
               res_rda$w %*% res_rda$q)
})

test_that("tt_fast_dl", {
  expect_equal(plsc(rmat1_s, rmat2_s)$d,
               tt_fast_dl(rmat1_s, rmat2_s, 'plsc')$d)
  expect_equal(cca(rmat1_s, rmat2_s)$d,
               tt_fast_dl(rmat1_s, rmat2_s, 'cca')$d)
  expect_equal(rda(rmat1_s, rmat2_s)$d,
               tt_fast_dl(rmat1_s, rmat2_s, 'rda')$d)
})

# library(bench)
#
# bench::mark(
# plsc(X_s, Y_s)$d,
# tt_fast_dl(X_s, Y_s, 'plsc', tol = .Machine$double.eps)$d
# )
#
# bench::mark(
#   cca(X_s, Y_s)$d,
#   tt_fast_dl(X_s, Y_s, 'cca', tol = .Machine$double.eps)$d
# )
#
# bench::mark(
#   rda(X_s, Y_s)$d,
#   tt_fast_dl(X_s, Y_s, 'rda', tol = .Machine$double.eps)$d
# )
