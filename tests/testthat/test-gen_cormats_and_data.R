ncolX <- 10
ncolY <- 5
nrows <- 100
seed_num <- 42

cor_mat <- multi_block_cor(cor_blocks = 0, ncolX = ncolX, ncolY = ncolY,
                           ncolX_blocks = 0, ncolY_blocks = 0,
                           cor_blocksX = 0, cor_blocksY = 0)

d_list <- gen_mvrnorm_from_cor_fast(cor_mat, nrows, ncolX, ncolY, seed = seed_num)
X <- d_list$X
Y <- d_list$Y

cor(X, Y)
