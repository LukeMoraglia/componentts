#' Compute pseudo-\eqn{R^2} values based on eigenvalues
#'
#' Given a vector of eigenvalues, this function computes the ratio between each
#'  eigenvalue and the sum of the eigenvalues after it. There are two types of
#'  pseudo-\eqn{R^2}s included: pseudo-\eqn{R^2} 1 includes only later eigenvalues
#'  in the sum, and pseudo-\eqn{R^2} 2 also includes the present eigenvalue in the
#'  sum.
#'
#' @param eigs a vector of eigenvalues for the computation
#' @param method "l+1" (default) corresponds to pseudo-\eqn{R^2} 1, and uses
#' all eigenvalues after the "current" eigenvalue as the divisor of the ratio.
#' "l" corresponds to pseudo-\eqn{R^2} 2, and uses all eigenvalues after the
#' "current" eigenvalue plus the current eigenvalue as the divisor of the ratio.
#'
#' @returns a vector the same length as `eigs` containing the pseudo-\eqn{R^2}s
#' @export
pseudo_R2 <- function(eigs, method = "l+1"){
  if(length(eigs) < 1){
    warning("pseudo_R2: eigs was length 0. returning `c(0)`")
    return(c(0))
  }

  if(method == "l+1"){
    divisors <- cumsum(eigs[length(eigs):1])[(length(eigs)-1):1]
    pseudo_R2s <- eigs / c(divisors, NA)
  }
  else if(method == "l"){
    divisors <- cumsum(eigs[length(eigs):1])[length(eigs):1]
    pseudo_R2s <- eigs / divisors
  }
  else{
    stop("pseudo_R2: unrecognized method")
  }

  return(pseudo_R2s)
}


#' Compute difference or ratio between consecutive values in a vector
#'
#' This function returns a vector the same length as its input vector `vec`
#' where each position has the following position subtracted from it or is
#' divided by it. This can be used with a vector of eigenvalues to get the
#' eigenvalue differences or eigenvalue ratios.
#'
#' @param vec an input vector of numerical values (such as eigenvalues)
#' @param method "diff" (default) or "ratio"
#'
#' @returns a vector the same length as `vec`. The last value is `NA`.
#' @export
diff_ratio <- function(vec, method = "diff"){
  if(length(vec) < 2){
    return(NA)
  }

  if(method == "diff"){
    res_vec <- c(vec[1:(length(vec)-1)] - vec[2:length(vec)], NA)
  }
  else if(method == "ratio"){
    res_vec <- c(vec[1:(length(vec)-1)] / vec[2:length(vec)], NA)
  }
  else{
    stop("diff_ratio: unrecognized method")
  }

  return(res_vec)
}

#' Compute RVDIM or RLSDIM from Dray (2008).
#'
#' This function computes for each eigenvalue an RVDIM value or an RLSDIM value
#' as described in Dray (2008).
#'
#' @inheritParams pseudo_R2
#' @param method "rvdim" (default) or "rlsdim"
#'
#' @returns a vector the same length as `eigs`
#' @export
rvdim <- function(eigs, method = "rvdim"){
  if(length(eigs) < 1){
    warning("rvdim: eigs was length 0. returning `c(0)`")
    return(c(0))
  }

  if(method == "rvdim"){
    divisors <- sqrt(cumsum((eigs^2)[length(eigs):1]))[length(eigs):1]
    res_vec <- eigs / divisors
  }
  else if(method == "rlsdim"){
    divisors <- sqrt(cumsum((eigs)[length(eigs):1]))[length(eigs):1]
    res_vec <- sqrt(eigs) / divisors
  }
  else{
    stop("rvdim: unrecognized method")
  }

  return(res_vec)

}

# Bartlett's statistic for each dimension As explained in Takane and Hwang (2002)
# Included only in case it's needed for back compatibility. Not maintained or tested.
bartlett <- function(eigs){
  cumsum((-log(1 - eigs))[length(eigs):1])[length(eigs):1]
}


