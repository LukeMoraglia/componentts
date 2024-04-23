#' Partial Least Squares Correlation
#'
#' A basic function for computing results of Partial Least Squares Correlation (PLSC).
#' Relies on the underlying [GSVD::gplssvd()] method.
#'
#' @param X first data matrix
#' @param Y second data matrix
#' @param tol (default: `.Machine$double.eps`) tolerance value for eigenvalues
#'
#' @returns
#' a list of results including:
#'
#'  * `d`: vector of singular values
#'  * `u`: matrix of left singular vectors
#'  * `v`: matrix of right singular vectors
#'  * `l`: vector of eigenvalues
#'  * `lx`: matrix of latent variables for `X`
#'  * `ly`: matrix of latent variables for `Y`
#'  * `p`: matrix of left generalized singular vectors. For PLSC, `p` and `u` are identical.
#'  * `q`: matrix of right generalized singular vectors. For PLSC, `q` and `v` are identical.
#'  * `m`: left metric matrix for the GSVD. For PLSC it is the identity matrix.
#'  * `w`: right metric matrix for the GSVD. For PLSC it is the identity matrix.
#'  * `f`: matrix of loadings for `X`. `lx` is computed by multiplying `X` and `f`. For PLSC,
#'         `f` and `p` are identical.
#'  * `g`: matrix of loadings for `Y`. `ly` is computed by multiplying `Y` and `g`. For PLSC,
#'         `g` and `q` are identical.
#'
#' @export
plsc <- function(X, Y, tol = .Machine$double.eps){

  res <- GSVD::gplssvd(X, Y, tol = tol)
  res$m <- diag(1, ncol(X), ncol(X))
  res$w <- diag(1, ncol(Y), ncol(Y))
  res$f <- res$p
  res$g <- res$q
  res$d_full <- NULL
  res$l_full <- NULL
  res$fi <- NULL
  res$fj <- NULL

  return(res)
}

#' Canonical Correlation Analysis
#'
#' A basic function for computing results of Canonical Correlation Analysis (CCA).
#' Relies on the underlying [GSVD::gplssvd()] method.
#'
#' @inheritParams plsc
#'
#' @returns a list of results including
#'  * `d`: vector of singular values
#'  * `u`: matrix of left singular vectors
#'  * `v`: matrix of right singular vectors
#'  * `l`: vector of eigenvalues
#'  * `lx`: matrix of latent variables for X
#'  * `ly`: matrix of latent variables for Y
#'  * `p`: matrix of left generalized singular vectors
#'  * `q`: matrix of right generalized singular vectors
#'  * `m`: left metric matrix for the GSVD. For CCA, it is the inverse crossproduct of `X`.
#'  * `w`: right metric matrix for the GSVD. For CCA, it is the inverse crossproduct of `Y`.
#'  * `f`: matrix of loadings for `X`. `lx` is computed by multiplying `X` and `f`. For CCA,
#'         `f` is `m` multiplied by `p`.
#'  * `g`: matrix of loadings for `Y`. `ly` is computed by multiplying `Y` and `g`. For CCA,
#'         `g` is `w` multiplied by `q`.
#'
#' @export
cca <- function(X, Y, tol = .Machine$double.eps){

  m <- MASS::ginv(crossprod(X))
  colnames(m) <- colnames(X)
  rownames(m) <- colnames(X)
  w <- MASS::ginv(crossprod(Y))
  colnames(w) <- colnames(Y)
  rownames(w) <- colnames(Y)
  res <- GSVD::gplssvd(X, Y,
                       XRW = m,
                       YRW = w,
                       tol = tol)
  res$m <- m
  res$w <- w
  res$f <- m %*% res$p
  res$g <- w %*% res$q
  res$d_full <- NULL
  res$l_full <- NULL
  res$fi <- NULL
  res$fj <- NULL

  return(res)
}

#' Redundancy Analysis
#'
#' A basic function for computing results of Redundancy Analysis (RDA).
#' Relies on the underlying [GSVD::gplssvd()] method.
#'
#' @inheritParams plsc
#'
#' @returns a list of results including
#'  * `d`: vector of singular values
#'  * `u`: matrix of left singular vectors
#'  * `v`: matrix of right singular vectors
#'  * `l`: vector of eigenvalues
#'  * `lx`: matrix of latent variables for X
#'  * `ly`: matrix of latent variables for Y
#'  * `p`: matrix of left generalized singular vectors
#'  * `q`: matrix of right generalized singular vectors
#'  * `m`: left metric matrix for the GSVD. For RDA, it is the inverse crossproduct of `X`.
#'  * `w`: right metric matrix for the GSVD. For RDA, it is the identity matrix.
#'  * `f`: matrix of loadings for `X`. `lx` is computed by multiplying `X` and `f`. For RDA,
#'         `f` is `m` multiplied by `p`.
#'  * `g`: matrix of loadings for `Y`. `ly` is computed by multiplying `Y` and `g`. For RDA,
#'         `g` and `q` are identical.
#'
#' @export
rda <- function(X, Y, tol = .Machine$double.eps){

  m <- MASS::ginv(crossprod(X))
  colnames(m) <- colnames(X)
  rownames(m) <- colnames(X)
  res <- GSVD::gplssvd(X, Y,
                       XRW = m,
                       tol = tol)
  res$m <- m
  res$w <- diag(1, ncol(Y), ncol(Y))
  colnames(res$w) <- colnames(Y)
  rownames(res$w) <- colnames(Y)
  res$f <- m %*% res$p
  res$g <- res$q
  res$d_full <- NULL
  res$l_full <- NULL
  res$fi <- NULL
  res$fj <- NULL

  return(res)
}

#' Select and run a two-table analysis: PLSC, CCA, or RDA
#'
#' This is a wrapper around [plsc()], [cca()], or [rda()].
#' The specific two-table analysis to run can be selected using `ttmethod`.
#'
#' @inheritParams plsc
#' @param ttmethod string denoting the two-table method to use,
#'  either 'plsc', 'cca', or 'rda'.
#'
#' @returns results of [plsc()], [cca()], or [rda()].
select_tt <- function(X, Y, ttmethod, tol = .Machine$double.eps){
  if(ttmethod == "plsc"){
    res <- plsc(X, Y, tol)
  }
  else if(ttmethod == "cca"){
    res <- cca(X, Y, tol)
  }
  else if(ttmethod == "rda"){
    res <- rda(X, Y, tol)
  }
  else{
    stop("select_tt: ttmethod is invalid")
  }
  return(res)
}

#' Fast singular values and eigenvalues for two-table methods
#'
#' Rather than compute all parts of the two-table methods,
#' `tt_fast_dl()` only computes the the singular values and eigenvalues.
#' Helpful for permutation tests when the goal is to get these values
#' with minimum computations.
#'
#' @inheritParams select_tt
#' @param Rx_invsqrt NULL (default) or
#'  a precomputed Rx inverse sqrt matrix (if you want to go faster)
#' @param Ry_invsqrt NULL (default) or
#'  a precomputed Ry inverse sqrt matrix (if you want to go faster)
#'
#' @returns a list with singular values (`d`) and eigenvalues (`l`)
#' @export
tt_fast_dl <- function(X, Y, ttmethod, tol = .Machine$double.eps, Rx_invsqrt = NULL, Ry_invsqrt = NULL){
  if(ttmethod == "plsc"){
    res_svd <- svd(mat_mult(t(X), Y), nu = 0, nv = 0)
  }
  else if(ttmethod == "cca"){
    if(is.null(Rx_invsqrt)){
      Rx_invsqrt <- fast_invsqrt_psd_matrix(mat_mult(t(X), X), symmetric = TRUE)
    }
    if(is.null(Ry_invsqrt)){
      Ry_invsqrt <- fast_invsqrt_psd_matrix(mat_mult(t(Y), Y), symmetric = TRUE)
    }

    matrices <- list(Rx_invsqrt, t(X), Y, Ry_invsqrt)
    res_svd <- svd(matrix_chain_multiplication(matrices), nu = 0, nv = 0)
  }
  else if(ttmethod == "rda"){
    if(is.null(Rx_invsqrt)){
      Rx_invsqrt <- fast_invsqrt_psd_matrix(mat_mult(t(X), X), symmetric = TRUE)
    }
    matrices <- list(Rx_invsqrt, t(X), Y)
    res_svd <- svd(matrix_chain_multiplication(matrices), nu = 0, nv = 0)
  }
  else{
    stop("tt_fast_dl: bad ttmethod")
  }
  res_svd$l <- res_svd$d^2
  svs_to_keep <- which(!(res_svd$l < tol))

  res_svd$d <- res_svd$d[svs_to_keep]
  res_svd$l <- res_svd$l[svs_to_keep]
  return(res_svd)
}
