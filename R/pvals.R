#' Compute p value given a statistic and a null distribution for the statistic
#'
#' This function computes p values based on null distributions from
#' permutation tests. It was written to be used with the results of
#' [SRP()], [DRP()], or [split_half()]. This function should not be used with results
#' from [multiple_hold()]; see [pvals_for_MH()] instead.
#'
#' @param stat a scalar or vector of observed statistics to test
#' @param stat_dist a matrix of size `length(stat)` by `n_iter` (number of permutations).
#'  Also accepts `n_iter` by `length(stat)`. For each value in `stat` there
#'  should be a corresponding row/column of `stat_dist` containing its null distribution.
#'  or `length(stat)` by `n_iter`.
#' @param test_type either "one_sided_right" (default), "one_sided_left", or "two_sided".
#'  This is the type of hypothesis test to do, aka one-tailed positive, one-tailed negative,
#'   or two-tailed. Essentially all two-table stopping rules use `one_sided_right`.
#' @param tol (default: 1e-6) tolerance used for comparisons between `stat` and `stat_dist`
#'
#' @returns a scalar or vector the same length as `stat` containing p values for `stat`
#'
#' @export
pval_from_dist <- function(stat, stat_dist, test_type = "one_sided_right", tol = 1e-6){
  stat <- as.vector(stat)
  stat_dist <- as.matrix(stat_dist)

  if(ncol(stat_dist) == nrow(stat_dist)){
    warning("pval_from_dist: stat_dist is square.
            Make sure that each column of stat_dist is a null distribution
            for the corresponding value in stat.")
  }
  # Make it so nrow(stat_dist) = length(stat)
  if(ncol(stat_dist) == length(stat)){
    stat_dist <- t(stat_dist)
  }
  else if(nrow(stat_dist) != length(stat)){
    stop("pval_from_dist: Either nrow(stat_dist) or ncol(stat_dist)
         must match length(stat)")
  }

  if(tolower(test_type) == "one_sided_right"){
    sig_count <- rowSums(stat < stat_dist + tol) # less than or equal to
  }
  else if(tolower(test_type) == "one_sided_left"){
    sig_count <- rowSums(stat + tol > stat_dist) # greater than or equal to
  }
  else if(tolower(test_type) == "two_sided"){
    sig_count <- rowSums(abs(stat) < abs(stat_dist) + tol) # greater than or equal to
  }
  else{
    stop("pval_from_dist: test_type unrecognized")
  }

  pval <- sig_count / ncol(stat_dist)
  return(pval)

}

#' Return p values for all statistics and null distributions in a list
#'
#' `pvals_all()` is a wrapper around [pval_from_dist()] that allows you to easily get
#' p values for several statistics and null distributions in a list, such as the
#' `perm_res` output of [SRP()], [DRP()], or [split_half()].
#'
#' @param dist_list list containing scalars/vectors that are statistics and
#'  matrices of null distributions. Every statistic must have an analagously named
#'  null distribution, so if there is a vector named "stat" there must be a matrix
#'  named "stat_dist" also in the list.
#' @inheritParams pval_from_dist
#'
#' @returns a list of length `length(dist_list)/2` containing p values for each
#' statistic/distribution in `dist_list`
#'
#' @export
pvals_all <- function(dist_list, test_type = "one_sided_right", tol = 1e-6){
  list_names <- names(dist_list)
  #should be even length
  if(length(list_names) %% 2 == 1){
    stop("pvals_all: dist_list should be even length
         (one item for each stat and stat distribution)")
  }

  stats_names <- list_names[!endsWith(list_names, "dist")]
  if(length(stats_names) != length(list_names)/2){
    stop("pvals_all: incorrect number of stats or distributions")
  }

  pvals_list <- list()

  for(stat in stats_names){
    pvals <- pval_from_dist(dist_list[[stat]], dist_list[[paste0(stat, "_dist")]],
                            test_type, tol)
    pvals_list[[paste0(stat, "_pvals")]] <- pvals
  }

  return(pvals_list)

}

#' Compute p values for the multiple holdout (MH) stopping rule
#'
#' This function takes the `perm_res` from the results of [multiple_hold()]
#' and uses it to compute p values in two different ways: the omnibus method
#' and the averaging method.
#'
#' @param dist_list the `perm_res` returned by [multiple_hold()] containing
#'  `rho` and `rho_dist`.
#' @param tol (default: 1e-6) tolerance for comparison between `rho` and `rho_dist`
#'
#' @returns list containing:
#'  * `pvals_omni`: a matrix of size `n_holdouts` by `length(rho)`.
#'  * `pvals_avg`: a vector the same length as `rho`
#'
#' @export
pvals_for_MH <- function(dist_list, tol = 1e-6){
  # Omnibus method
  rho = dist_list$rho
  rho_dist = aperm(dist_list$rho_dist, c(3, 2, 1))
  tf_brick <- abind::abind(apply(rho_dist, MARGIN = 3,
                    FUN = function(rho_sheet){rho < rho_sheet + tol}, # less than or equal to
                    simplify = FALSE), along = 3)
  sig_count <- apply(tf_brick, c(1,2), sum)
  pvals_omni <- sig_count / dim(rho_dist)[3]


  # Averaging method
  rho_means <- colMeans(rho)
  rho_means_dist <- apply(rho_dist, 3, colMeans)
  sig_count <- rowSums(rho_means < rho_means_dist + tol) # less than or equal to
  pvals_avg <- sig_count / ncol(rho_means_dist)

  return(list(pvals_omni = pvals_omni,
              pvals_avg = pvals_avg))
}
