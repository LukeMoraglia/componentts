#' Inverse square root of positive semi-definite matrix
#'
#' This function was stolen from the GSVD package by Derek Beaton
#' and modified to eliminate the input checks
#' and allow for additional arguments to be passed to `tolerance_eigen` and by extension
#' `eigen`. This makes it much faster by passing the `symmetric` argument when
#' matrix symmetry is known.
#'
#' @param x positive semi-definite matrix
#' @param ... everything else passed to `tolerance_eigen`
#'
#' @return the inverse square root matrix
#' @export
fast_invsqrt_psd_matrix <- function(x, ...){
  ## tolerance_eigen
  res <- GSVD::tolerance_eigen(x, tol = 1e-13, ...)

  ## rebuild
  return(t(t(res$vectors) * (1/sqrt(res$values)) ) %*% t(res$vectors))

}


#' Center and normalize the columns of a matrix to a sum of squares of 1.
#'
#' Centers each column and scales each column to have a sum of squares of 1.
#'
#' @param X matrix to center and scale
#'
#' @return centered and scaled X
#' @export
scale_SS1 <- function(X){
  centering <- apply(X, 2, mean)
  sqrt_ss <- function(col){return(sqrt(sum((col - mean(col))^2)))}
  scaling <- apply(X, 2, sqrt_ss)
  X_scaled <- scale(X, center = centering, scale = scaling)
  return(X_scaled)
}
