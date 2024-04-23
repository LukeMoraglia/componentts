set.seed(42)
X <- matrix(sample.int(20, 20*3, replace = TRUE), 20, 3)
set.seed(69)
Y <- matrix(sample.int(20, 20*4, replace = TRUE), 20, 4)
X_s <- scale_SS1(X)
Y_s <- scale_SS1(Y)

set.seed(24)
srp_plsc1 <- SRP(X_s, Y_s, 'plsc', rand_type = 'rand_rows', n_iter = 399, parallel = FALSE)

res1 <- pval_from_dist(srp_plsc1$perm_res$l_prop, srp_plsc1$perm_res$l_prop_dist)
res2 <- pvals_all(srp_plsc1$perm_res)

res_mh_plsc <- multiple_hold(X_s, Y_s, "plsc", percent_train = 0.6)

res3 <- pvals_for_MH(res_mh_plsc$perm_res)
