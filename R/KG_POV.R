#' Kaiser-Guttman Criterion (Kaiser Line)
#'
#' A simple stopping rule to keep components with an eigenvalue greater than the average
#' eigenvalue. This average eigenvalue is the Kaiser line or `k_line`. The function
#' can either take `X`, `Y`, and `ttmethod` and compute the results, or a vector
#' of already computed eigenvalues `eigs` can be passed in. This stopping rule is
#' not recommended.
#'
#' @inheritParams select_tt
#' @param X (default: `NULL`) first data matrix
#' @param Y (default: `NULL`) second data matrix
#' @param ttmethod (default: `NULL`) string denoting the two-table method to use,
#'  either 'plsc', 'cca', or 'rda'.
#' @param eigs (default: `NULL`) a vector of eigenvalues for the test.
#'
#' @returns list containing:
#'  * `l`: vector of eigenvalues
#'  * `k_line`: the average of the eigenvalues, above this value is considered
#'    important or significant.
#'  * `imp_comp`: boolean vector saying whether each component is above the
#'    `k_line` or not.
#'
#' @export
kaiser_guttman <- function(X = NULL, Y = NULL, ttmethod = NULL,
                           eigs = NULL, tol = .Machine$double.eps){
  if(!is.null(eigs)){ # if eigs is provided, ignore X, Y, and ttmethod
    l <- eigs
  }
  else if(!is.null(X) & !is.null(Y) & !is.null(ttmethod)){ # X, Y, and ttmethod provided
    fixed_res <- select_tt(X, Y, ttmethod, tol)
    l <- fixed_res$l
  }
  else{
    stop("kaiser_guttman: you must provide eigs or provide X, Y, and ttmethod")
  }

  k_line <- mean(l)
  imp_comp <- l > k_line

  return(list(l = l,
              k_line = k_line,
              imp_comp = imp_comp))
}


#' Percentage of variance (POV) explained
#'
#' A simple stopping rule which returns the cumulative percentage of variance
#' or 'inertia' explained by the components of a two-table analysis. One can then
#' use the cumulative percentage to keep the first set of components which reach
#' a certain threshold, such as 80%. The function can either take `X`, `Y`, and
#'`ttmethod` and compute the results, or a vector of already computed
#' eigenvalues `eigs` can be passed in. This stopping rule is not recommended.
#'
#' @inheritParams kaiser_guttman
#'
#' @returns list containing:
#'  * `l`: the vector of eigenvalues from the analysis
#'  * `inertia`: the total inertia (variance), which is the sum of the eigenvalues
#'  * `perc_inertia`: a vector the same length as `l` giving the percentage of
#'    inertia explained by each component.
#'  * `cum_perc_inertia`: a vector the same length as `l` giving the cumulative
#'    percentage of inertia explained from the first component up to that component
#'
#' @export
pov <- function(X = NULL, Y = NULL, ttmethod = NULL,
                eigs = NULL, tol = .Machine$double.eps){
  if(!is.null(eigs)){ # if eigs is provided, ignore X, Y, and ttmethod
    l <- eigs
  }
  else if(!is.null(X) & !is.null(Y) & !is.null(ttmethod)){ # X, Y, and ttmethod provided
    fixed_res <- select_tt(X, Y, ttmethod, tol)
    l <- fixed_res$l
  }
  else{
    stop("pov: you must provide eigs or provide X, Y, and ttmethod")
  }

  inertia <- sum(l)
  perc_inertia <- l / inertia
  cum_perc_inertia <- cumsum(perc_inertia)

  return(list(l = l,
              inertia = inertia,
              perc_inertia = perc_inertia,
              cum_perc_inertia = cum_perc_inertia))
}
