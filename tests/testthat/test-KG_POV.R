set.seed(42)
X <- matrix(sample.int(20, 20*3, replace = TRUE), 20, 3)
set.seed(69)
Y <- matrix(sample.int(20, 20*4, replace = TRUE), 20, 4)
X_s <- scale_SS1(X)
Y_s <- scale_SS1(Y)

res_KG <- kaiser_guttman(X_s, Y_s, 'plsc')
res_pov <- pov(X_s, Y_s, 'plsc')
