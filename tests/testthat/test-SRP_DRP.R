# Just want to test that these run without errors

set.seed(42)
X <- matrix(sample.int(20, 20*3, replace = TRUE), 20, 3)
set.seed(69)
Y <- matrix(sample.int(20, 20*4, replace = TRUE), 20, 4)
X_s <- scale_SS1(X)
Y_s <- scale_SS1(Y)

set.seed(24)
#srp_plsc1_par <- SRP(X_s, Y_s, 'plsc', rand_type = 'rand_rows', n_iter = 399, parallel = TRUE)
set.seed(24)
srp_plsc1 <- SRP(X_s, Y_s, 'plsc', rand_type = 'rand_rows', n_iter = 399, parallel = FALSE)

srp_plsc2 <- SRP(X_s, Y_s, 'plsc', rand_type = 'rand_each_col')

srp_cca1 <- SRP(X_s, Y_s, 'cca', rand_type = 'rand_rows')
srp_cca2 <- SRP(X_s, Y_s, 'cca', rand_type = 'rand_each_col')

srp_rda1 <- SRP(X_s, Y_s, 'rda', rand_type = 'rand_rows')
srp_rda2 <- SRP(X_s, Y_s, 'rda', rand_type = 'rand_each_col')

set.seed(24)
drp_plsc1 <- DRP(X_s, Y_s, 'plsc', rand_type = 'rand_rows')
drp_plsc2 <- DRP(X_s, Y_s, 'plsc', rand_type = 'rand_each_col')

drp_cca1 <- DRP(X_s, Y_s, 'cca', rand_type = 'rand_rows')
drp_cca2 <- DRP(X_s, Y_s, 'cca', rand_type = 'rand_each_col')

drp_rda1 <- DRP(X_s, Y_s, 'rda', rand_type = 'rand_rows')
drp_rda2 <- DRP(X_s, Y_s, 'rda', rand_type = 'rand_each_col')




