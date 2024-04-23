# The functions in this file are wrappers around their respective stopping rule
# functions. They get the results, time how long it takes, and compute the pvals.
# Everything is returned as a list of results.
run_SRP <- function(X, Y, ttmethod, n_true_comp = NULL, seed_num = NULL, cor_mat_params = NULL,
                    n_iter = 99, rand_type = "rand_rows",
                    tol = .Machine$double.eps,
                    parallel = FALSE){
  set.seed(seed_num)
  start_time <- Sys.time()
  stop_rule_res <- SRP(X, Y, ttmethod, n_iter, rand_type, tol, parallel)
  end_time <- Sys.time()
  pvals <- pvals_all(stop_rule_res$perm_res)
  # store the metadata on this run
  run_res <- list(stop_rule_res = stop_rule_res,
                  pvals = pvals,
                  cor_mat_params = cor_mat_params,
                  seed_num = seed_num,
                  ttmethod = ttmethod,
                  n_true_comp = n_true_comp,
                  n_iter = n_iter,
                  rand_type = rand_type,
                  tol = tol,
                  runtime = end_time - start_time)
  return(run_res)
}

run_DRP <- function(X, Y, ttmethod, n_true_comp = NULL, seed_num = NULL,
                    cor_mat_params = NULL, n_iter = 99, rand_type = "rand_rows",
                    tol = .Machine$double.eps, parallel = FALSE,
                    n_dim_comp = NULL){
  set.seed(seed_num)
  start_time <- Sys.time()
  stop_rule_res <- DRP(X, Y, ttmethod, n_iter, rand_type, tol, parallel,
                       n_dim_comp)
  end_time <- Sys.time()
  pvals <- pvals_all(stop_rule_res$perm_res)
  # store the metadata on this run
  run_res <- list(stop_rule_res = stop_rule_res,
                  pvals = pvals,
                  cor_mat_params = cor_mat_params,
                  seed_num = seed_num,
                  ttmethod = ttmethod,
                  n_true_comp = n_true_comp,
                  n_iter = n_iter,
                  rand_type = rand_type,
                  tol = tol,
                  runtime = end_time - start_time)
  return(run_res)
}

run_MH <- function(X, Y, ttmethod, n_true_comp = NULL, seed_num = NULL,
                   cor_mat_params = NULL, n_iter = 99, rand_type = "rand_rows",
                   tol = .Machine$double.eps,
                   n_holdouts = 10, percent_train = 0.9, parallel = FALSE,
                   n_dim_comp = NULL){
  set.seed(seed_num)
  start_time <- Sys.time()
  stop_rule_res <- multiple_hold(X, Y, ttmethod, n_iter, rand_type,
                                 tol, n_holdouts, percent_train, parallel,
                                 n_dim_comp)
  end_time <- Sys.time()
  # need a custom pval computing method for MH
  pvals <- pvals_for_MH(stop_rule_res$perm_res)

  # store the metadata on this run
  run_res <- list(stop_rule_res = stop_rule_res,
                  pvals = pvals,
                  cor_mat_params = cor_mat_params,
                  seed_num = seed_num,
                  ttmethod = ttmethod,
                  n_true_comp = n_true_comp,
                  n_iter = n_iter,
                  rand_type = rand_type,
                  tol = tol,
                  runtime = end_time - start_time)
  return(run_res)
}

run_SHR <- function(X, Y, ttmethod, n_true_comp = NULL, seed_num = NULL, cor_mat_params = NULL,
                    n_iter = 99, rand_type = "rand_rows",
                    tol = .Machine$double.eps, n_splits = 200){
  set.seed(seed_num)
  start_time <- Sys.time()
  stop_rule_res <- split_half(X, Y, ttmethod, n_iter, rand_type, tol, n_splits)
  end_time <- Sys.time()
  pvals <- pvals_all(stop_rule_res$perm_res)
  # store the metadata on this run
  run_res <- list(stop_rule_res = stop_rule_res,
                  pvals = pvals,
                  cor_mat_params = cor_mat_params,
                  seed_num = seed_num,
                  ttmethod = ttmethod,
                  n_true_comp = n_true_comp,
                  n_iter = n_iter,
                  rand_type = rand_type,
                  tol = tol,
                  runtime = end_time - start_time)
  return(run_res)
}
