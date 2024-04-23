#' Multiple holdout (MH)
#'
#' MH is a two-table stopping rule that splits data into a training and holdout
#' set, assesses how the holdout data does on results of the training data,
#' and compares that to permuted data using a permutation test with deflation.
#' The function does not compute p values, but that can be accomplished by passing the
#' results to [pvals_for_MH()].
#'
#' @inheritParams DRP
#' @param n_holdouts (default: 10) the number of holdout sets to compute
#' @param percent_train (default: .9) the percentage of observations used
#'  for the training set
#'
#' @returns a list containing:
#'  * `fixed_res`: a list containing results from the original, unpermuted analysis.
#'  The results of [select_tt()] are passed through unchanged.
#'  * `perm_res`: a list containing results from the permutation iterations.
#'    * `rho`: a matrix of size `n_holdouts` by number of components (`n_dim`), containing
#'    the rho values. Rho assesses the agreement between the latent variables computed
#'    from the holdout set, so there are `n_dim` rho values for each holdout set.
#'    * `rho_dist`: an array of size `n_iter + 1` by `n_dim` by `n_holdouts` containing
#'    the null distributions of the rho values.
#'
#' @export
multiple_hold <- function(X, Y, ttmethod, n_iter = 99, rand_type = 'rand_rows',
                          tol = .Machine$double.eps, n_holdouts = 10,
                          percent_train = .9, parallel = FALSE, n_dim_comp = NULL){
  fixed_res <- select_tt(X, Y, ttmethod, tol)
  n_train <- round(nrow(X) * percent_train)
  n_hold <- nrow(X) - n_train

  if(n_hold < 3){
    warning("multiple_hold: holdout set only has 2 rows, changing to 3 rows. This will change percent_train slightly.")
    n_hold <- 3
    n_train <- nrow(X) - n_hold
  }

  n_dim <- min(ncol(X), ncol(Y), n_train)

  if(is.null(n_dim_comp)){
    n_dim <- min(ncol(X), ncol(Y), n_train)
  }
  else if(n_dim_comp > 0 & n_dim_comp <= n_dim){
    n_dim <- n_dim_comp
  }
  else{
    stop("multiple_hold: n_dim_comp argument invalid. should be integer or NULL")
  }


  rho <- matrix(NA, nrow = n_holdouts, ncol = n_dim)
  rho_dist <- array(NA, dim = c(n_iter + 1, n_dim, n_holdouts))

  if(parallel){
    n_cores <- parallel::detectCores() - 1
    clus <- parallel::makeCluster(n_cores,
                                  type = "PSOCK")
    doParallel::registerDoParallel(clus)
  }

  if(tolower(rand_type) == 'rr'){
    rand_type <- 'rand_rows'
  }
  else if(tolower(rand_type) == 'rec'){
    rand_type <- 'rand_each_col'
  }

  for(dim in 1:n_dim){
    def <- deflate_XY(X, Y, ttmethod, fixed_res, dim)

    X_def <- def$X_def
    Y_def <- def$Y_def

    for(holdout in 1:n_holdouts){
      ind_train <- sample(1:nrow(X), n_train, replace = FALSE)
      ind_hold <- c(1:nrow(X))[!(1:nrow(X) %in% ind_train)]
      X_train <- X_def[ind_train,]
      Y_train <- Y_def[ind_train,]
      X_hold <- X_def[ind_hold,]
      Y_hold <- Y_def[ind_hold,]


      # precompute Rx_invsqrt and Ry_invsqrt once if rand_type is rand_rows,
      # since they stay the same after permutation
      Rx_invsqrt <- NULL
      Ry_invsqrt <- NULL
      if(tolower(rand_type) == 'rand_rows'){
        if(ttmethod == 'cca'){
          Rx_invsqrt <- fast_invsqrt_psd_matrix(mat_mult(t(X_train), X_train), symmetric = TRUE)
          Ry_invsqrt <- fast_invsqrt_psd_matrix(mat_mult(t(Y_train), Y_train), symmetric = TRUE)
        }
        else if(ttmethod == 'rda'){
          Rx_invsqrt <- fast_invsqrt_psd_matrix(mat_mult(t(X_train), X_train), symmetric = TRUE)
        }
      }

      rho_hold <- mh_fast_rho(X_train, Y_train, X_hold, Y_hold, ttmethod, tol, Rx_invsqrt, Ry_invsqrt)
      rho[holdout, dim] <- rho_hold

      if(parallel){
        # main permutation loop
        rand_res_list <- foreach::foreach(i = 1:n_iter) %dorng% {
          if(tolower(rand_type) == 'rand_rows'){
            X_rand <- rand_rows_cpp(X_train)
            Y_rand <- Y_train
          }
          else if(tolower(rand_type) == 'rand_each_col'){
            X_rand <- rand_each_col_cpp(X_train)
            Y_rand <- rand_each_col_cpp(Y_train)
          }
          else{
            stop("multiple_hold: rand_type is invalid")
          }

          mh_fast_rho(X_rand, Y_rand, X_hold, Y_hold, ttmethod, tol, Rx_invsqrt, Ry_invsqrt)
        }
      }
      else{ # not parallel
        # main permutation loop
        rand_res_list <- foreach::foreach(i = 1:n_iter) %do% {
          if(tolower(rand_type) == 'rand_rows'){
            X_rand <- rand_rows_cpp(X_train)
            Y_rand <- Y_train
          }
          else if(tolower(rand_type) == 'rand_each_col'){
            X_rand <- rand_each_col_cpp(X_train)
            Y_rand <- rand_each_col_cpp(Y_train)
          }
          else{
            stop("multiple_hold: rand_type is invalid")
          }

          mh_fast_rho(X_rand, Y_rand, X_hold, Y_hold, ttmethod, tol, Rx_invsqrt, Ry_invsqrt)
        }
      }
      for(i in 1:n_iter){
        rho_dist[i, dim, holdout] <- rand_res_list[[i]]
      }

      rho_dist[n_iter + 1, dim, holdout] <- rho_hold
    }
  }

  if(parallel){
    parallel::stopCluster(clus)
  }

  perm_res <- list(rho = rho,
                   rho_dist = rho_dist)
  return(list(fixed_res = fixed_res,
              perm_res = perm_res))
}

# Helper function for `multiple_hold()` to quickly compute rho
mh_fast_rho <- function(X_train, Y_train, X_hold, Y_hold, ttmethod,
                        tol = .Machine$double.eps, Rx_invsqrt = NULL, Ry_invsqrt = NULL){
  if(ttmethod == "plsc"){
    res_svd <- svd(mat_mult(t(X_train), Y_train), nu = 1, nv = 1)
    lvx1 <- mat_mult(X_hold, res_svd$u)
    lvy1 <- mat_mult(Y_hold, res_svd$v)
  }
  else if(ttmethod == "cca"){
    if(is.null(Rx_invsqrt)){
      Rx_invsqrt <- fast_invsqrt_psd_matrix(mat_mult(t(X_train), X_train), symmetric = TRUE)
    }
    if(is.null(Ry_invsqrt)){
      Ry_invsqrt <- fast_invsqrt_psd_matrix(mat_mult(t(Y_train), Y_train), symmetric = TRUE)
    }
    matrices <- list(Rx_invsqrt, t(X_train), Y_train, Ry_invsqrt)
    res_svd <- svd(matrix_chain_multiplication(matrices), nu = 1, nv = 1)
    lvx1 <- mat_mult(X_hold, mat_mult(Rx_invsqrt, res_svd$u))
    lvy1 <- mat_mult(Y_hold, mat_mult(Ry_invsqrt, res_svd$v))
  }
  else if(ttmethod == "rda"){
    if(is.null(Rx_invsqrt)){
      Rx_invsqrt <- fast_invsqrt_psd_matrix(mat_mult(t(X_train), X_train), symmetric = TRUE)
    }
    matrices <- list(Rx_invsqrt, t(X_train), Y_train)
    res_svd <- svd(matrix_chain_multiplication(matrices), nu = 1, nv = 1)
    lvx1 <- mat_mult(X_hold, mat_mult(Rx_invsqrt, res_svd$u))
    lvy1 <- mat_mult(Y_hold, res_svd$v)
  }
  else{
    stop("mh_fast_rho: bad ttmethod")
  }

  return(abs(stats::cor(as.vector(lvx1), as.vector(lvy1))))

}

#' Split half reliability (SHR)
#'
#' SHR is a two-table stopping rule that splits data in half, assesses how well each
#' half does on its own, and compares that to permuted data using a permutation test
#' without deflation.
#' The function does not actually compute the p values for each component; that is accomplished by
#' passing the results to [pval_from_dist()] or [pvals_all()].
#'
#'
#' @inheritParams SRP
#' @param n_splits (default: 100) number of data splits used to compute the `Pcorr_mean` and `Qcorr_mean`.
#'
#' @returns a list containing:
#'  * `fixed_res`: a list containing results from the original, unpermuted analysis.
#'    * `d`: singular values
#'    * `u`: left singular vectors
#'    * `v`: right singular vectors
#'    * `f`: loadings for `X`
#'    * `g`: loadings for `Y`
#'  * `perm_res`: a list containing results from the permutation iterations.
#'    * `Pcorr_mean`: a vector containing the mean of the `Pcorr` values for each component.
#'    * `Pcorr_mean_dist`: a matrix of size `n_iter + 1` by `length(Pcorr_mean)` containing the
#'    null distributions of the Pcorr values
#'    * `Qcorr_mean`: a vector containing the mean of the `Qcorr` values for each component.
#'    * `Qcorr_mean_dist`: a matrix of size `n_iter + 1` by `length(Qcorr_mean)` containing the
#'    null distributions of the Qcorr values
#'
#' @export
split_half <- function(X, Y, ttmethod, n_iter = 99, rand_type = 'rand_rows',
                       tol = .Machine$double.eps, n_splits = 100){
  fixed_res <- shr_fast_res(X, Y, ttmethod, tol)
  n_dim <- length(fixed_res$d)
  n_per_split <- round(nrow(X) / 2)

  d_inv <- diag(fixed_res$d^(-1))
  g_d_inv <- mat_mult(fixed_res$g, d_inv)
  f_d_inv <- mat_mult(fixed_res$f, d_inv)

  Pcorr <- matrix(NA, n_splits, n_dim)
  Qcorr <- matrix(NA, n_splits, n_dim)

  for(split in 1:n_splits){
    ind1 <- sample(1:nrow(X), n_per_split, replace = FALSE)
    ind2 <- c(1:nrow(X))[!(1:nrow(X) %in% ind1)]
    X1 <- X[ind1,]
    Y1 <- Y[ind1,]
    X2 <- X[ind2,]
    Y2 <- Y[ind2,]
    R1 <- mat_mult(t(X1), Y1)
    R2 <- mat_mult(t(X2), Y2)

    P1 <- mat_mult(R1, g_d_inv)
    P2 <- mat_mult(R2, g_d_inv)

    Q1 <- mat_mult(t(R1), f_d_inv)
    Q2 <- mat_mult(t(R2), f_d_inv)

    Pcorr[split,] <- abs(diag(stats::cor(P1, P2)))
    Qcorr[split,] <- abs(diag(stats::cor(Q1, Q2)))
  }

  Pcorr_mean <- apply(Pcorr, 2, mean)
  Qcorr_mean <- apply(Qcorr, 2, mean)

  Pcorr_mean_dist <- matrix(NA, n_iter + 1, n_dim)
  Qcorr_mean_dist <- matrix(NA, n_iter + 1, n_dim)

  if(tolower(rand_type) == 'rr'){
    rand_type <- 'rand_rows'
  }
  else if(tolower(rand_type) == 'rec'){
    rand_type <- 'rand_each_col'
  }

  Rx_invsqrt <- NULL
  Ry_invsqrt <- NULL
  if(tolower(rand_type) == 'rand_rows'){
    if(ttmethod == 'cca'){
      Rx_invsqrt <- fast_invsqrt_psd_matrix(mat_mult(t(X), X), symmetric = TRUE)
      Ry_invsqrt <- fast_invsqrt_psd_matrix(mat_mult(t(Y), Y), symmetric = TRUE)
    }
    else if(ttmethod == 'rda'){
      Rx_invsqrt <- fast_invsqrt_psd_matrix(mat_mult(t(X), X), symmetric = TRUE)
    }
  }

  for(i in 1:n_iter){
    if(tolower(rand_type) == 'rand_rows'){
      X_rand <- rand_rows(X)
      Y_rand <- Y
    }
    else if(tolower(rand_type) == 'rand_each_col'){
      X_rand <- rand_each_col(X)
      Y_rand <- rand_each_col(Y)
    }
    else{
      stop("split_half: rand_type is invalid")
    }

    rand_res <- shr_fast_res(X_rand, Y_rand, ttmethod, tol, Rx_invsqrt, Ry_invsqrt)

    if(length(rand_res$d) < n_dim){ # rand_res dim is too small
      rand_res$d[(length(rand_res$d)+1):n_dim] <- 0
    }
    else if(length(rand_res$d) > n_dim){
      rand_res$d <- rand_res$d[1:n_dim]
    }
    d_inv_rand <- diag(rand_res$d^(-1))

    ind1 <- sample(1:nrow(X_rand), n_per_split, replace = FALSE)
    ind2 <- c(1:nrow(X_rand))[!(1:nrow(X_rand) %in% ind1)]
    X1 <- X_rand[ind1,]
    Y1 <- Y_rand[ind1,]
    X2 <- X_rand[ind2,]
    Y2 <- Y_rand[ind2,]
    R1 <- mat_mult(t(X1), Y1)
    R2 <- mat_mult(t(X2), Y2)
    g_d_inv <- mat_mult(rand_res$g, d_inv_rand)
    f_d_inv <- mat_mult(rand_res$f, d_inv_rand)
    P1 <- mat_mult(R1, g_d_inv)
    P2 <- mat_mult(R2, g_d_inv)

    Q1 <- mat_mult(t(R1), f_d_inv)
    Q2 <- mat_mult(t(R2), f_d_inv)

    Pcorr_mean_dist[i,] <- abs(diag(stats::cor(P1, P2)))
    Qcorr_mean_dist[i,] <- abs(diag(stats::cor(Q1, Q2)))
  }

  Pcorr_mean_dist[n_iter + 1,] <- Pcorr_mean
  Qcorr_mean_dist[n_iter + 1,] <- Qcorr_mean

  perm_res <- list(Pcorr_mean = Pcorr_mean,
                   Pcorr_mean_dist = Pcorr_mean_dist,
                   Qcorr_mean = Qcorr_mean,
                   Qcorr_mean_dist = Qcorr_mean_dist)
  return(list(fixed_res = fixed_res,
              perm_res = perm_res))
}

# Helper function to quickly compute the values required for SHR
shr_fast_res <- function(X, Y, ttmethod, tol = .Machine$double.eps,
                         Rx_invsqrt = NULL, Ry_invsqrt = NULL){
  if(ttmethod == "plsc"){
    res_svd <- svd(mat_mult(t(X), Y))
    svs_to_keep <- which(!(res_svd$d^2 < tol))
    res_svd$d <- res_svd$d[svs_to_keep]
    res_svd$f <- res_svd$u[,svs_to_keep]
    res_svd$g <- res_svd$v[,svs_to_keep]
  }
  else if(ttmethod == "cca"){
    if(is.null(Rx_invsqrt)){
      Rx_invsqrt <- fast_invsqrt_psd_matrix(mat_mult(t(X), X), symmetric = TRUE)
    }
    if(is.null(Ry_invsqrt)){
      Ry_invsqrt <- fast_invsqrt_psd_matrix(mat_mult(t(Y), Y), symmetric = TRUE)
    }

    matrices <- list(Rx_invsqrt, t(X), Y, Ry_invsqrt)

    res_svd <- svd(matrix_chain_multiplication(matrices))
    svs_to_keep <- which(!(res_svd$d^2 < tol))
    res_svd$d <- res_svd$d[svs_to_keep]
    res_svd$f <- mat_mult(Rx_invsqrt, res_svd$u[,svs_to_keep])
    res_svd$g <- mat_mult(Ry_invsqrt, res_svd$v[,svs_to_keep])
  }
  else if(ttmethod == "rda"){
    if(is.null(Rx_invsqrt)){
      Rx_invsqrt <- fast_invsqrt_psd_matrix(mat_mult(t(X), X), symmetric = TRUE)
    }
    matrices <- list(Rx_invsqrt, t(X), Y)
    res_svd <- svd(matrix_chain_multiplication(matrices))
    svs_to_keep <- which(!(res_svd$d^2 < tol))
    res_svd$d <- res_svd$d[svs_to_keep]
    res_svd$f <- mat_mult(Rx_invsqrt, res_svd$u[,svs_to_keep])
    res_svd$g <- res_svd$v[,svs_to_keep]
  }
  else{
    stop("shr_fast_res: bad ttmethod")
  }

  return(res_svd)

}
