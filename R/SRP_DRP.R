#' Simultaneous randomization procedure (SRP)
#'
#' SRP is a two-table stopping rule that runs a simple permutation test without deflation.
#' The function computes several test statistics and their null distributions.
#' It does not actually compute the p values for each component; that is accomplished by
#' passing the results to [pval_from_dist()] or [pvals_all()].
#'
#' @inheritParams select_tt
#' @param n_iter (default: 99) number of permutation iterations
#' @param rand_type (default: 'rand_rows') how to randomize the matrices,
#'    either 'rand_rows'/'RR' or 'rand_each_col'/'REC'.
#'    The RR option permutes entire rows of the `X` matrix while leaving `Y` unchanged.
#'    The REC option permutes every column of both `X` and `Y`.
#' @param parallel (default: FALSE) use a parallelized loop to do the main permutations.
#'    This will use the `parallel` package to utilize multiple cores simultaneously
#'    (using [parallel::detectCores()] - 1 cores).
#'
#' @import doRNG
#' @import foreach
#'
#' @returns a list containing:
#'  * `fixed_res`: a list containing results from the original, unpermuted analysis
#'    * `d`: singular values from the two-table analysis
#'    * `l`: eigenvalues from the two-table analysis
#'  * `perm_res`: a list containing results from the permutation iterations.
#'  For each test statistic, there is a vector of the observed test statistics
#'  for each component, and there is a matrix of the null distribution for each
#'  component obtained from the permutation iterations.
#'    * `d`: a vector of singular values
#'    * `d_dist`: a matrix of size `n_iter + 1` by `length(d)` containing the
#'    null distributions of the singular values
#'    * `l_prop`: a vector of the eigenvalue proportions (EVP). EVP is computed as
#'    a component's eigenvalue divided by the sum of all the eigenvalues.
#'    * `l_prop_dist`: a matrix of size `n_iter + 1` by `length(l_prop)` containing the
#'    null distributions of the EVPs.
#'    * `pseudo_R2_1`: a vector of pseudo-\eqn{R^2} 1 values. Pseudo-\eqn{R^2} 1 is
#'    computed as a component's eigenvalue divided by the sum of all the following
#'    eigenvalues, in contrast to pseudo-\eqn{R^2} 2 below.
#'    The last component is assigned `NA`.
#'    * `pseudo_R2_1_dist`: a matrix of size `n_iter + 1` by `length(pseudo_R2_1)`
#'    containing the null distributions of the pseudo-\eqn{R^2} 1 values. The last
#'    column corresponding to the last component are all `NA`.
#'    * `pseudo_R2_2`: a vector of pseudo-\eqn{R^2} 2 values. Pseudo-\eqn{R^2} 2 is
#'    computed as a component's eigenvalue divided by the sum of itself and all the following
#'    eigenvalues, in contrast to pseudo-\eqn{R^2} 1 above.
#'    The last component is assigned `NA`.
#'    * `pseudo_R2_2_dist`: a matrix of size `n_iter + 1` by `length(pseudo_R2_2)`
#'    containing the null distributions of the pseudo-\eqn{R^2} 2 values. The last
#'    column corresponding to the last component are all `NA`.
#'    * `l_diff`: a vector containing the eigenvalue differences, the difference between
#'    a component's eigenvalue and the following eigenvalue. The last component is assigned `NA`.
#'    * `l_diff_dist`: a matrix of size `n_iter + 1` by `length(l_diff)` containing the
#'    null distributions of the eigenvalue differences. The last
#'    column corresponding to the last component are all `NA`.
#'    * `l_ratio`: a vector containing the eigenvalue ratios, the ratio between
#'    a component's eigenvalue and the following eigenvalue. The last component is assigned `NA`.
#'    * `l_ratio_dist`: a matrix of size `n_iter + 1` by `length(l_ratio)` containing the
#'    null distributions of the eigenvalue ratios. The last
#'    column corresponding to the last component are all `NA`.
#'    * `rvdim`: a vector containing the RVDIM values, a value proposed by Dray (2008) that
#'    assesses the similarity between the rank-one matrix created by a component, and
#'    the matrix created by all of the following components. The last component is assigned `NA`.
#'    * `rvdim_dist`: a matrix of size `n_iter + 1` by `length(rvdim)` containing the
#'    null distributions of the RVDIM values. The last
#'    column corresponding to the last component are all `NA`.
#'
#' @export
SRP <- function(X, Y, ttmethod, n_iter = 99, rand_type = 'rand_rows',
                tol = .Machine$double.eps, parallel = FALSE){
  fixed_res <- tt_fast_dl(X, Y, ttmethod, tol)
  n_dim <- length(fixed_res$d)
  d_dist <- matrix(NA, nrow = n_iter + 1, ncol = n_dim)
  l_prop_dist <- matrix(NA, n_iter + 1, n_dim)
  pseudo_R2_1_dist <- matrix(NA, n_iter + 1, n_dim)
  pseudo_R2_2_dist <- matrix(NA, n_iter + 1, n_dim)
  l_diff_dist <- matrix(NA, n_iter + 1, n_dim)
  l_ratio_dist <- matrix(NA, n_iter + 1, n_dim)
  rvdim_dist <- matrix(NA, n_iter + 1, n_dim)

  if(tolower(rand_type) == 'rr'){
    rand_type <- 'rand_rows'
  }
  else if(tolower(rand_type) == 'rec'){
    rand_type <- 'rand_each_col'
  }

  # precompute Rx_invsqrt and Ry_invsqrt once if rand_type is rand_rows, since
  # they stay the same after permutation
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

  if(parallel){
    n_cores <- parallel::detectCores() - 1
    clus <- parallel::makeCluster(n_cores,
                                  type = "PSOCK")
    doParallel::registerDoParallel(clus)

    # main permutation loop
    rand_res_list <- foreach::foreach(i = 1:n_iter) %dorng% {
      if(tolower(rand_type) == 'rand_rows'){
        X_rand <- rand_rows_cpp(X)
        Y_rand <- Y
      }
      else if(tolower(rand_type) == 'rand_each_col'){
        X_rand <- rand_each_col_cpp(X)
        Y_rand <- rand_each_col_cpp(Y)
      }
      else{
        stop("SRP: rand_type is invalid")
      }

      tt_fast_dl(X_rand, Y_rand, ttmethod, tol, Rx_invsqrt, Ry_invsqrt)
    }

    parallel::stopCluster(clus)
  }
  else{
    # main permutation loop
    rand_res_list <- foreach::foreach(i = 1:n_iter) %do% {
      if(tolower(rand_type) == 'rand_rows'){
        X_rand <- rand_rows_cpp(X)
        Y_rand <- Y
      }
      else if(tolower(rand_type) == 'rand_each_col'){
        X_rand <- rand_each_col_cpp(X)
        Y_rand <- rand_each_col_cpp(Y)
      }
      else{
        stop("SRP: rand_type is invalid")
      }

      tt_fast_dl(X_rand, Y_rand, ttmethod, tol, Rx_invsqrt, Ry_invsqrt)
    }
  }

  for(i in 1:n_iter){
    rand_res <- rand_res_list[[i]]
    if(length(rand_res$d) < n_dim){ # rand_res dim is too small
      rand_res$d[(length(rand_res$d)+1):n_dim] <- 0
      rand_res$l[(length(rand_res$l)+1):n_dim] <- 0
    }
    d_dist[i, ] <- rand_res$d[1:n_dim]
    l_prop_dist[i, ] <- rand_res$l[1:n_dim] / sum(rand_res$l[1:n_dim])
    pseudo_R2_1_dist[i, ] <- pseudo_R2(rand_res$l[1:n_dim], method = "l+1")
    pseudo_R2_2_dist[i, ] <- pseudo_R2(rand_res$l[1:n_dim], method = "l")
    l_diff_dist[i, ] <- diff_ratio(rand_res$l[1:n_dim], method = "diff")
    l_ratio_dist[i, ] <- diff_ratio(rand_res$l[1:n_dim], method = "ratio")
    rvdim_dist[i, ] <- rvdim(rand_res$l[1:n_dim], method = "rvdim")
  }

  d <- fixed_res$d
  d_dist[n_iter + 1,] <- d
  l_prop <- fixed_res$l / sum(fixed_res$l)
  l_prop_dist[n_iter + 1, ] <- l_prop
  pseudo_R2_1 <- pseudo_R2(fixed_res$l, method = "l+1")
  pseudo_R2_1_dist[n_iter + 1, ] <- pseudo_R2_1
  pseudo_R2_1_dist[, n_dim] <- NA
  pseudo_R2_2 <- pseudo_R2(fixed_res$l, method = "l")
  pseudo_R2_2_dist[n_iter + 1, ] <- pseudo_R2_2
  pseudo_R2_2_dist[, n_dim] <- NA
  l_diff <- diff_ratio(fixed_res$l, method = "diff")
  l_diff_dist[n_iter + 1, ] <- l_diff
  l_diff_dist[, n_dim] <- NA
  l_ratio <- diff_ratio(fixed_res$l, method = "ratio")
  l_ratio_dist[n_iter + 1, ] <- l_ratio
  l_ratio_dist[, n_dim] <- NA
  rvdim <- rvdim(fixed_res$l, method = "rvdim")
  rvdim_dist[n_iter + 1, ] <- rvdim
  rvdim_dist[, n_dim] <- NA

  perm_res <- list(d = d,
                   d_dist = d_dist,
                   l_prop = l_prop,
                   l_prop_dist = l_prop_dist,
                   pseudo_R2_1 = pseudo_R2_1,
                   pseudo_R2_1_dist = pseudo_R2_1_dist,
                   pseudo_R2_2 = pseudo_R2_2,
                   pseudo_R2_2_dist = pseudo_R2_2_dist,
                   l_diff = l_diff,
                   l_diff_dist = l_diff_dist,
                   l_ratio = l_ratio,
                   l_ratio_dist = l_ratio_dist,
                   rvdim = rvdim,
                   rvdim_dist = rvdim_dist)


  return(list(fixed_res = fixed_res,
              perm_res = perm_res)
  )
}



#' Deflation randomization procedure (DRP)
#'
#' DRP is a two-table stopping rule that runs a sequential permutation test with deflation.
#' Each component is tested individually, and the effects of previous components are
#' removed using deflation.
#' The function computes several test statistics and their null distributions.
#' It does not actually compute the p values for each component; that is accomplished by
#' passing the results to [pval_from_dist()] or [pvals_all()].
#'
#' @inheritParams SRP
#' @param n_dim_comp (default: `NULL`) the number of dimensions to compute null distributions for.
#'    `NULL` means all dimensions are tested.
#'    Testing a smaller number of components will cut down on computational time.
#'
#' @import doRNG
#' @import foreach
#'
#' @returns See [SRP()] for the list of test statistics and null distributions that are returned,
#'  with the only difference that `DRP()` does not include `l_prop` or `l_prop_dist`.
#'
#' @export
DRP <- function(X, Y, ttmethod, n_iter = 99, rand_type = 'rand_rows',
                tol = .Machine$double.eps, parallel = FALSE, n_dim_comp = NULL){
  fixed_res <- select_tt(X, Y, ttmethod, tol)
  n_dim <- length(fixed_res$d)

  if(is.null(n_dim_comp)){
    n_dim <- length(fixed_res$d)
  }
  else if(n_dim_comp > 0 & n_dim_comp <= n_dim){
    n_dim <- n_dim_comp
  }
  else{
    stop("DRP: n_dim_comp argument invalid. should be integer or NULL")
  }

  d_dist <- matrix(NA, nrow = n_iter + 1, ncol = n_dim)
  pseudo_R2_1_dist <- matrix(NA, n_iter + 1, n_dim)
  pseudo_R2_2_dist <- matrix(NA, n_iter + 1, n_dim)
  l_diff_dist <- matrix(NA, n_iter + 1, n_dim)
  l_ratio_dist <- matrix(NA, n_iter + 1, n_dim)
  rvdim_dist <- matrix(NA, n_iter + 1, n_dim)

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

    # precompute Rx_invsqrt and Ry_invsqrt once if rand_type is rand_rows, since
    # they stay the same after permutation
    Rx_invsqrt <- NULL
    Ry_invsqrt <- NULL
    if(tolower(rand_type) == 'rand_rows'){
      if(ttmethod == 'cca'){
        Rx_invsqrt <- fast_invsqrt_psd_matrix(mat_mult(t(X_def), X_def), symmetric = TRUE)
        Ry_invsqrt <- fast_invsqrt_psd_matrix(mat_mult(t(Y_def), Y_def), symmetric = TRUE)
      }
      else if(ttmethod == 'rda'){
        Rx_invsqrt <- fast_invsqrt_psd_matrix(mat_mult(t(X_def), X_def), symmetric = TRUE)
      }
    }

    if(parallel){
      # main permutation loop
      rand_res_list <- foreach::foreach(i = 1:n_iter) %dorng% {
        if(tolower(rand_type) == 'rand_rows'){
          X_rand <- rand_rows_cpp(X_def)
          Y_rand <- Y_def
        }
        else if(tolower(rand_type) == 'rand_each_col'){
          X_rand <- rand_each_col_cpp(X_def)
          Y_rand <- rand_each_col_cpp(Y_def)
        }
        else{
          stop("DRP: rand_type is invalid")
        }

        tt_fast_dl(X_rand, Y_rand, ttmethod, tol, Rx_invsqrt, Ry_invsqrt)
      }
    }
    else{ # not parallel
      # main permutation loop
      rand_res_list <- foreach::foreach(i = 1:n_iter) %do% {
        if(tolower(rand_type) == 'rand_rows'){
          X_rand <- rand_rows_cpp(X_def)
          Y_rand <- Y_def
        }
        else if(tolower(rand_type) == 'rand_each_col'){
          X_rand <- rand_each_col_cpp(X_def)
          Y_rand <- rand_each_col_cpp(Y_def)
        }
        else{
          stop("DRP: rand_type is invalid")
        }

        tt_fast_dl(X_rand, Y_rand, ttmethod, tol, Rx_invsqrt, Ry_invsqrt)
      }
    }

    for(i in 1:n_iter){
      rand_res <- rand_res_list[[i]]
      d_dist[i, dim] <- rand_res$d[1]
      pseudo_R2_1_dist[i, dim] <- pseudo_R2(rand_res$l, method = "l+1")[1]
      pseudo_R2_2_dist[i, dim] <- pseudo_R2(rand_res$l, method = "l")[1]
      l_diff_dist[i, dim] <- diff_ratio(rand_res$l, method = "diff")[1]
      l_ratio_dist[i, dim] <- diff_ratio(rand_res$l, method = "ratio")[1]
      rvdim_dist[i, dim] <- rvdim(rand_res$l, method = "rvdim")[1]
    }
  }


  if(parallel){
    parallel::stopCluster(clus)
  }

  d <- fixed_res$d[1:n_dim]
  d_dist[n_iter + 1,] <- d
  pseudo_R2_1 <- pseudo_R2(fixed_res$l, method = "l+1")[1:n_dim]
  pseudo_R2_1_dist[n_iter + 1, ] <- pseudo_R2_1
  pseudo_R2_2 <- pseudo_R2(fixed_res$l, method = "l")[1:n_dim]
  pseudo_R2_2_dist[n_iter + 1, ] <- pseudo_R2_2
  l_diff <- diff_ratio(fixed_res$l, method = "diff")[1:n_dim]
  l_diff_dist[n_iter + 1, ] <- l_diff
  l_ratio <- diff_ratio(fixed_res$l, method = "ratio")[1:n_dim]
  l_ratio_dist[n_iter + 1, ] <- l_ratio
  rvdim <- rvdim(fixed_res$l, method = "rvdim")[1:n_dim]
  rvdim_dist[n_iter + 1, ] <- rvdim


  if(n_dim == length(fixed_res$d)){ # when n_dim is the last possible dim, these are NA
    pseudo_R2_1_dist[, n_dim] <- NA
    pseudo_R2_2_dist[, n_dim] <- NA
    l_diff_dist[, n_dim] <- NA
    l_ratio_dist[, n_dim] <- NA
    rvdim_dist[, n_dim] <- NA
  }

  perm_res <- list(d = d,
                   d_dist = d_dist,
                   pseudo_R2_1 = pseudo_R2_1,
                   pseudo_R2_1_dist = pseudo_R2_1_dist,
                   pseudo_R2_2 = pseudo_R2_2,
                   pseudo_R2_2_dist = pseudo_R2_2_dist,
                   l_diff = l_diff,
                   l_diff_dist = l_diff_dist,
                   l_ratio = l_ratio,
                   l_ratio_dist = l_ratio_dist,
                   rvdim = rvdim,
                   rvdim_dist = rvdim_dist)

  return(list(fixed_res = fixed_res,
              perm_res = perm_res)
  )
}

# Deflate X and Y of the effects of earlier dimensions.
# This function is a helper function for DRP and MH.
# Given the results `res` and a dimension `dim`, it removes the
# effects of all dimensions/components before `dim`.
deflate_XY <- function(X, Y, ttmethod, res, dim){
  if(ttmethod == "plsc"){
    wX <- as.matrix(res$u[,0:(dim-1)])
    projX <- diag(1, ncol(X), ncol(X)) - (wX %*% t(wX))
    X_def <- X %*% projX

    wY <- as.matrix(res$v[,0:(dim-1)])
    projY <- diag(1, ncol(Y), ncol(Y)) - (wY %*% t(wY))
    Y_def <- Y %*% projY
  }
  else if(ttmethod == "cca"){
    R_x <- t(X) %*% X
    wX <- as.matrix((GSVD::invsqrt_psd_matrix(R_x) %*% res$u)[,0:(dim-1)])
    projX <- diag(1, ncol(X), ncol(X)) - (wX %*% t(wX) %*% R_x)
    X_def <- X %*% projX

    R_y <- t(Y) %*% Y
    wY <- as.matrix((GSVD::invsqrt_psd_matrix(R_y) %*% res$v)[,0:(dim-1)])
    projY <- diag(1, ncol(Y), ncol(Y)) - (wY %*% t(wY) %*% R_y)
    Y_def <- Y %*% projY
  }
  else if(ttmethod == "rda"){
    R_x <- t(X) %*% X
    wX <- as.matrix((GSVD::invsqrt_psd_matrix(R_x) %*% res$u)[,0:(dim-1)])
    projX <- diag(1, ncol(X), ncol(X)) - (wX %*% t(wX) %*% R_x)
    X_def <- X %*% projX

    wY <- as.matrix(res$v[,0:(dim-1)])
    projY <- diag(1, ncol(Y), ncol(Y)) - (wY %*% t(wY))
    Y_def <- Y %*% projY
  }

  return(list(X_def = X_def, Y_def = Y_def))
}


