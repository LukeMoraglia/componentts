#' Create a correlation matrix between X and Y with multiple blocks of variables
#'
#' This function creates a joint correlation matrix for the variables of X and Y,
#' with the main feature being the blocks of correlated variables placed on the
#' diagonal.
#'
#' @param cor_blocks a vector of correlation values corresponding to each block.
#'  `length(cor_blocks)` determines the number of blocks (`n_blocks`).
#' @param ncolX integer number of columns that belong to the X matrix
#' @param ncolY integer number of columns that belong to the Y matrix
#' @param ncolX_blocks vector same length as `cor_blocks` that gives the
#'  number of X columns in each block.
#' @param ncolY_blocks vector same length as `cor_blocks` that gives the
#'  number of Y columns in each block.
#' @param cor_blocksX vector the same length as `cor_blocks` that contains the
#'  within X correlations for each block
#' @param cor_blocksY vector the same length as `cor_blocks` that contains the
#'  within X correlations for each block
#' @param offXY scalar correlation between variables of X and Y that aren't in a block
#' @param offX scalar correlation between variables of X that aren't in a block
#' @param offY scalar correlation between variables of Y that aren't in a block
#'
#' @returns a correlation matrix based on the parameters, with size
#'  `ncolX + ncolY` by `ncolX + ncolY`. The first `ncolX` rows/columns
#'  correspond to the X variables, and the last `ncolY` rows/columns to the
#'  Y variables.
#'
#' @export
multi_block_cor <- function(cor_blocks,
                            ncolX, ncolY,
                            ncolX_blocks, ncolY_blocks,
                            cor_blocksX, cor_blocksY,
                            offXY = 0, offX = 0, offY = 0){
  ncolXY <- ncolX + ncolY
  n_blocks <- length(cor_blocks)
  cor_mat <- matrix(offXY, ncolXY, ncolXY)


  # X start ind
  X_s_ind <- 1
  # Y start index
  Y_s_ind <- ncolX + 1
  cor_mat[1:ncolX, 1:ncolX] <- offX
  cor_mat[Y_s_ind:ncolXY, Y_s_ind:ncolXY] <- offY

  for(b in 1:n_blocks){
    X_e_ind <- X_s_ind + ncolX_blocks[b] - 1
    Y_e_ind <- Y_s_ind + ncolY_blocks[b] - 1
    cor_mat[X_s_ind:X_e_ind, Y_s_ind:Y_e_ind] <- cor_blocks[b]
    cor_mat[X_s_ind:X_e_ind, X_s_ind:X_e_ind] <- cor_blocksX[b]
    cor_mat[Y_s_ind:Y_e_ind, Y_s_ind:Y_e_ind] <- cor_blocksY[b]
    X_s_ind <- X_e_ind + 1
    Y_s_ind <- Y_e_ind + 1
  }
  diag(cor_mat) <- 1
  cor_mat[lower.tri(cor_mat)] <- t(cor_mat)[lower.tri(cor_mat)]
  return(cor_mat)
}

#' Generate multivariate normal data from a correlation matrix.
#'
#' This function samples two tables of data, `X` and `Y`,
#' from a multivariate normal distribution using the joint correlation matrix
#' between `X` and `Y`, likely created using [multi_block_cor()].
#'
#' @param cor_mat a matrix of size `ncolX + ncolY` by `ncolX + ncolY`. This is
#'  the joint correlation matrix of `X` and `Y`, likely created using [multi_block_cor()].
#' @param nrows integer number of rows to sample
#' @param ncolX integer number of columns for `X`. Rows/columns 1:ncolX of `cor_mat`
#'  correspond with `X`.
#' @param ncolY integer number of columns for `X`. Rows/columns (ncolX+1):(ncolX+ncolY)
#'  of `cor_mat` correspond with `Y`.
#' @param mu vector of means for each column. default is `NULL` which sets all means
#'  to zero.
#' @param seed (default:NULL) optional seed for reproducibility. Because of the
#'  underlying C++ code, the normal R `set.seed()` will not work and it must be
#'  passed to the function instead.
#'
#' @returns list containing:
#'  * `X`: matrix of size `nrows` by `ncolX`
#'  * `Y`: matrix of size `nrows` by `ncolY`
#' @export
gen_mvrnorm_from_cor_fast <- function(cor_mat, nrows, ncolX, ncolY, mu = NULL, seed = NULL){
  ncolXY <- ncolX + ncolY
  if(is.null(mu)){
    mu <- rep(0, ncolXY)
  }
  data <- my_rmvnorm(nrows, mu, cor_mat, seed)
  X <- data[,1:ncolX]
  Y <- data[,(ncolX + 1):ncolXY]
  return(list(X = X, Y = Y))
}

# Rather than using Rfast::rmvnorm(), this is a custom version that
# runs faster.
my_rmvnorm <- function(n, mu, sigma, seed = NULL)
{
  p <- length(mu)
  if (!is.null(seed))
    RcppZiggurat::zsetseed(seed)
  x <- Rfast::matrnorm(n, p)
  Rfast::mat.mult(x, chol(sigma)) + rep(mu, rep(n, p))
}
