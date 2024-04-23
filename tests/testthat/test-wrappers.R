# To get accurate benchmarks with compiled code, do not use devtools::load_all()
# Restart R, reinstall package, restart R, and load package with library(componentts)
library(tictoc)

ncolX <- 100
ncolY <- 50
nrows <- 1000
seed_num <- 42

cor_mat <- multi_block_cor(cor_blocks = 0, ncolX = ncolX, ncolY = ncolY,
                           ncolX_blocks = 0, ncolY_blocks = 0,
                           cor_blocksX = 0, cor_blocksY = 0)
d_list <- gen_mvrnorm_from_cor_fast(cor_mat, nrows, ncolX, ncolY, seed = seed_num)
X <- d_list$X
Y <- d_list$Y

X_s <- scale_SS1(X)
Y_s <- scale_SS1(Y)

# Commented out to prevent a slow run of tests during check()

# bench::mark(
#   srp_res1 = componentts:::run_SRP(X_s, Y_s, "plsc", n_true_comp = 0,
#                                    seed_num = seed_num, cor_mat_params = NULL,
#                                    n_iter = 199, rand_type = "rand_rows", parallel = FALSE),
#   iterations = 4,
#   filter_gc = FALSE
# )
# # 7/19/23 1.6s
# # 7/19/23 671ms (after adding mat_mult and nu=nv=0)
# # 7/20/23 essentially the same after adding rand_rows_cpp
# # 8/17/23 adding parallel makes it go to 4.7s hahaha
#
# bench::mark(
#   srp_res2 = componentts:::run_SRP(X_s, Y_s, "cca", n_true_comp = 0,
#                                    seed_num = seed_num, cor_mat_params = NULL,
#                                    n_iter = 199, rand_type = "rand_rows"),
#   iterations = 4
# )
# # 7/19/23 667ms
# # 7/20/23 essentially the same after adding rand_rows_cpp
#
# bench::mark(
#   srp_res3 = componentts:::run_SRP(X_s, Y_s, "rda", n_true_comp = 0,
#                                    seed_num = seed_num, cor_mat_params = NULL,
#                                    n_iter = 199, rand_type = "rand_rows"),
#   iterations = 4
# )
# # 7/19/23 672ms
# # 7/20/23 essentially the same after adding rand_rows_cpp
#
# bench::mark(
#   srp_res4 = componentts:::run_SRP(X_s, Y_s, "plsc", n_true_comp = 0,
#                                    seed_num = seed_num, cor_mat_params = NULL,
#                                    n_iter = 199, rand_type = "rand_each_col"),
#   iterations = 3
# )
# # 7/19/23 4.85s
# # 7/19/23 3.5s (after adding mat_mult and nu=nv=0)
# # 7/20/23 1.3s (after rand_each_col_cpp)
#
# bench::mark(
#   srp_res5 = componentts:::run_SRP(X_s, Y_s, "cca", n_true_comp = 0,
#                                    seed_num = seed_num, cor_mat_params = NULL,
#                                    n_iter = 199, rand_type = "rand_each_col"),
#   iterations = 3
# )
# # 7/19/23 5.32s
# # 7/20/23 4.71s (after rand_each_col_cpp)
# # 7/20/23 3.9s (get rid of crossprod)
#
# bench::mark(
#   srp_res6 = componentts:::run_SRP(X_s, Y_s, "rda", n_true_comp = 0,
#                                    seed_num = seed_num, cor_mat_params = NULL,
#                                    n_iter = 199, rand_type = "rand_each_col"),
#   iterations = 3
# )
# # 7/19/23 5.25s
# # 7/20/23 4.1s (after rand_each_col_cpp)
# # 7/20/23 3.0s (get rid of crossprod)
#
# profvis::profvis({
#   componentts:::run_SRP(X_s, Y_s, "plsc", n_true_comp = 0,
#                         seed_num = seed_num, cor_mat_params = NULL,
#                         n_iter = 399, rand_type = "rand_rows", parallel = FALSE)
# })
#
# profvis::profvis({
#   componentts:::run_SRP(X_s, Y_s, "cca", n_true_comp = 0,
#                         seed_num = seed_num, cor_mat_params = NULL,
#                         n_iter = 399, rand_type = "rand_rows")
# })
#
# profvis::profvis({
#   componentts:::run_SRP(X_s, Y_s, "plsc", n_true_comp = 0,
#                         seed_num = seed_num, cor_mat_params = NULL,
#                         n_iter = 399, rand_type = "rand_each_col")
# })
#
# profvis::profvis({
#   componentts:::run_SRP(X_s, Y_s, "cca", n_true_comp = 0,
#                         seed_num = seed_num, cor_mat_params = NULL,
#                         n_iter = 999, rand_type = "rand_each_col", parallel = TRUE)
# })
#
#
# ### DRP
#
# bench::mark(
#   drp_res1 = componentts:::run_DRP(X_s, Y_s, "plsc", n_true_comp = 0,
#                                    seed_num = seed_num, cor_mat_params = NULL,
#                                    n_iter = 49, rand_type = "rand_rows", parallel = FALSE),
#   iterations = 1
# )
# # 8/15/23 9.4s
#
#
# tic()
# drp_res2 = componentts:::run_DRP(X_s, Y_s, "cca", n_true_comp = 0,
#                                  seed_num = seed_num, cor_mat_params = NULL,
#                                  n_iter = 49, rand_type = "rand_rows", parallel = FALSE)
# toc()
#
# # 8/15/23 11.9s
# # 9/27/23 13.7s (no parallel), 17.1s (with parallel)
#
# tic()
# drp_res2 = componentts:::run_DRP(X_s, Y_s, "cca", n_true_comp = 0,
#                                  seed_num = seed_num, cor_mat_params = NULL,
#                                  n_iter = 99, rand_type = "rand_rows", parallel = TRUE)
# toc()
# # 9/27/23 25.8s (no parallel) 22.16 (with parallel)
#
# bench::mark(
#   drp_res3 = componentts:::run_DRP(X_s, Y_s, "rda", n_true_comp = 0,
#                                    seed_num = seed_num, cor_mat_params = NULL,
#                                    n_iter = 49, rand_type = "rand_rows"),
#   iterations = 1
# )
#
# # 8/15/23 13.3s
#
# bench::mark(
#   drp_res4 = componentts:::run_DRP(X_s, Y_s, "plsc", n_true_comp = 0,
#                                    seed_num = seed_num, cor_mat_params = NULL,
#                                    n_iter = 49, rand_type = "rand_each_col"),
#   iterations = 1
# )
# # 8/15/23 16.4s
#
#
# tic()
# drp_res5 = componentts:::run_DRP(X_s, Y_s, "cca", n_true_comp = 0,
#                                  seed_num = seed_num, cor_mat_params = NULL,
#                                  n_iter = 49, rand_type = "rand_each_col", parallel = TRUE)
# toc()
# # 8/15/23 39.5s
# # 9/27/23 42.3s (no parallel), 23.5s (with parallel)!
#
#
# bench::mark(
#   drp_res6 = componentts:::run_DRP(X_s, Y_s, "rda", n_true_comp = 0,
#                                    seed_num = seed_num, cor_mat_params = NULL,
#                                    n_iter = 49, rand_type = "rand_each_col"),
#   iterations = 1
# )
#
# # 8/15/23 32.2s
#
# profvis::profvis({
#   componentts:::run_DRP(X_s, Y_s, "plsc", n_true_comp = 0,
#                         seed_num = seed_num, cor_mat_params = NULL,
#                         n_iter = 49, rand_type = "rand_rows")
# })
#
# profvis::profvis({
#   componentts:::run_DRP(X_s, Y_s, "cca", n_true_comp = 0,
#                         seed_num = seed_num, cor_mat_params = NULL,
#                         n_iter = 49, rand_type = "rand_rows")
# })
#
# profvis::profvis({
#   componentts:::run_DRP(X_s, Y_s, "plsc", n_true_comp = 0,
#                         seed_num = seed_num, cor_mat_params = NULL,
#                         n_iter = 49, rand_type = "rand_each_col")
# })
#
# profvis::profvis({
#   componentts:::run_DRP(X_s, Y_s, "cca", n_true_comp = 0,
#                         seed_num = seed_num, cor_mat_params = NULL,
#                         n_iter = 49, rand_type = "rand_each_col")
# })
#
#
# MH

# tic()
# mh_res1 = componentts:::run_MH(X_s, Y_s, "plsc", n_true_comp = 0,
#                                  seed_num = seed_num, cor_mat_params = NULL,
#                                  n_iter = 49, rand_type = "rand_rows", parallel = FALSE)
# toc()
# # 145s (no parallel), 114s (with parallel)
#
# bench::mark(
#   mh_res2 = componentts:::run_MH(X_s, Y_s, "cca", n_true_comp = 0,
#                                  seed_num = seed_num, cor_mat_params = NULL,
#                                  n_iter = 19, rand_type = "rand_rows", n_holdouts = 5),
#   iterations = 1
# )
# # 8/17/23 31.4s
#
# profvis::profvis({
#   componentts:::run_MH(X_s, Y_s, "plsc", n_true_comp = 0,
#                        seed_num = seed_num, cor_mat_params = NULL,
#                        n_iter = 19, rand_type = "rand_rows")
# })

# SHR

# tic()
# shr_res1 = componentts:::run_SHR(X_s, Y_s, 'plsc', n_true_comp = 0,
#                                  seed_num = seed_num, cor_mat_params = NULL,
#                                  n_iter = 199, rand_type = 'rand_rows', n_splits = 200)
# toc()
# # 9/27/23 4.22s, 3.6s (mat_mult), 2.9 (shr_fast optimized)
#
# tic()
# shr_res2 = componentts:::run_SHR(X_s, Y_s, 'cca', n_true_comp = 0,
#                                  seed_num = seed_num, cor_mat_params = NULL,
#                                  n_iter = 199, rand_type = 'rand_rows', n_splits = 200)
# toc()
# # 9/27/23 10.6s, 7.5s (mat_mult), 2.81 (shr_fast optimized)
#
# tic()
# shr_res3 = componentts:::run_SHR(X_s, Y_s, 'cca', n_true_comp = 0,
#                                  seed_num = seed_num, cor_mat_params = NULL,
#                                  n_iter = 199, rand_type = 'rand_each_col', n_splits = 200)
# toc()
# # 9/27/23 11.87s, 10.6 (mat_mult), 7.2 (shr_fast optimized)
#
# profvis::profvis({
#   shr_res2 = componentts:::run_SHR(X_s, Y_s, 'cca', n_true_comp = 0,
#                                    seed_num = seed_num, cor_mat_params = NULL,
#                                    n_iter = 199, rand_type = 'rand_rows', n_splits = 200)
# })
