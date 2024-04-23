#' Randomize (permute) the rows of a matrix
#'
#' This function randomizes or permutes a matrix by reordering complete rows of
#' the matrix. This results in a matrix where the relationships (correlations)
#' within the matrix are preserved. But when compared to another matrix with different
#' variables and the same participants, the relationships between the two matrices
#' are now broken.
#'
#' This function is the option used when selecting "rand_rows" in the stopping rule
#' functions such as [SRP()] or [DRP()]. However, all the stopping rule functions
#' in this package use a faster version written with `Rcpp`.
#'
#' @param X the matrix to permute
#'
#' @returns a matrix with the rows of `X` reordered
#' @export
rand_rows <- function(X){
  new_row_order <- sample(1:nrow(X), replace = FALSE)
  return(X[new_row_order,])
}

#' Randomize (permute) each column of a matrix
#'
#' This function randomizes or permutes a matrix by reordering the values in each
#' column of the matrix. This results in a matrix where the relationships (correlations)
#' within the matrix are broken, and when compared to another matrix with different
#' variables and the same participants, the relationships between the two matrices
#' are also now broken.
#'
#' This function is the option used when selecting "rand_each_col" in the stopping rule
#' functions such as [SRP()] or [DRP()]. However, all the stopping rule functions
#' in this package use a faster version written with `Rcpp`.
#'
#' @param X the matrix to permute
#'
#' @returns a matrix with the values in each column of `X` reordered
#' @export
rand_each_col <- function(X){
  return(apply(X, 2, sample, replace = FALSE))
}
