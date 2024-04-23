set.seed(42)
X <- matrix(sample.int(20, 20*3, replace = TRUE), 20, 3)
set.seed(69)
Y <- matrix(sample.int(20, 20*4, replace = TRUE), 20, 4)
X_s <- scale_SS1(X)
Y_s <- scale_SS1(Y)


res_plsc <- plsc(X_s, Y_s)
res_cca <- cca(X_s, Y_s)
res_rda <- rda(X_s, Y_s)

# lv_plsc <- lvs_from_res(X_s, Y_s, "plsc", res_plsc)
# lv_cca <- lvs_from_res(X_s, Y_s, "cca", res_cca)
# lv_rda <- lvs_from_res(X_s, Y_s, "rda", res_rda)
#
# test_that("latent variable computation", {
#   expect_equal(lv_plsc$lvx, res_plsc$lx)
#   expect_equal(lv_plsc$lvy, res_plsc$ly)
#   expect_equal(lv_cca$lvx, res_cca$lx)
#   expect_equal(lv_cca$lvy, res_cca$ly)
#   expect_equal(lv_rda$lvx, res_rda$lx)
#   expect_equal(lv_rda$lvy, res_rda$ly)
# })

#debugonce(multiple_hold)
res_mh_plsc <- multiple_hold(X_s, Y_s, "plsc", percent_train = 0.6)
res_mh_cca <- multiple_hold(X_s, Y_s, 'cca')
res_mh_rda <- multiple_hold(X_s, Y_s, 'rda')

res_sh_plsc <- split_half(X_s, Y_s, "plsc")
res_sh_cca <- split_half(X_s, Y_s, 'cca')
res_sh_rda <- split_half(X_s, Y_s, 'rda')

rho_plsc <- mh_fast_rho(X_s, Y_s, X_s, Y_s, "plsc")
rho_cca <- mh_fast_rho(X_s, Y_s, X_s, Y_s, "cca")
rho_rda <- mh_fast_rho(X_s, Y_s, X_s, Y_s, "rda")

test_that("mh_fast_rho", {
  expect_equal(rho_plsc,
               abs(cor(res_plsc$lx[,1], res_plsc$ly[,1])))
  expect_equal(rho_cca,
               abs(cor(res_cca$lx[,1], res_cca$ly[,1])))
  expect_equal(rho_rda,
               abs(cor(res_rda$lx[,1], res_rda$ly[,1])))
})

pvals_for_MH(res_mh_plsc$perm_res)
