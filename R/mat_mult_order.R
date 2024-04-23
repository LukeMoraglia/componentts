# Helper function to compute the cost and split matrices for the matrix
# chain multiplication algorithm. See for instance, the pseudocode given
# at https://en.wikipedia.org/wiki/Matrix_chain_multiplication.
matrix_chain_order <- function(p) {
  n <- length(p) - 1

  # Initialize cost and split matrices
  cost <- matrix(0, n, n)
  split <- matrix(0, n, n)

  # Compute cost and split matrices
  for (l in 2:n) {
    for (i in 1:(n-l+1)) {
      j <- i + l - 1
      cost[i, j] <- Inf
      for (k in i:(j-1)) {
        q <- cost[i, k] + cost[k+1, j] + p[i]*p[k+1]*p[j+1]
        if (q < cost[i, j]) {
          cost[i, j] <- q
          split[i, j] <- k
        }
      }
    }
  }

  # Construct optimal multiplication order
  optimal_order <- construct_optimal_order(split, 1, n)

  # Return cost, split, and optimal order
  return(list(cost = cost, split = split, optimal_order = optimal_order))
}

# Helper function to construct optimal multiplication order
construct_optimal_order <- function(split, i, j) {
  if (i == j) {
    return(paste0("A", i))
  } else {
    k <- split[i, j]
    left <- construct_optimal_order(split, i, k)
    right <- construct_optimal_order(split, k+1, j)
    return(paste0("(", left, "*", right, ")"))
  }
}

# Helper function to get the dimensions of the matrices and pass them to
# `matrix_chain_order`
get_split_mat <- function(matrices){
  n <- length(matrices)

  # Extract dimensions of matrices
  dimensions <- sapply(matrices, function(m) dim(m)[1])
  dimensions <- c(dimensions, dim(matrices[[n]])[2])

  # Find optimal multiplication order
  return(matrix_chain_order(dimensions)$split)
}


#' Matrix Chain Multiplication
#'
#' Given a list of compatible matrices to multiply, this function returns their
#' multiplication, using an efficient order of matrix multiplications. The order
#' of matrix multiplications is computed using a matrix chain multiplication
#' algorithm (see for instance: https://en.wikipedia.org/wiki/Matrix_chain_multiplication).
#'
#' @param matrices list of compatible matrices to multiply. Matrices are multiplied
#'  in the order they are given in the list.
#' @param split (default: NULL) the split matrix computed by `get_split_mat` can
#'  be supplied here, or it is computed automatically if `split = NULL`.
#'
#' @returns a matrix that is the result of multiplying the matrices in `matrices`
#' @export
matrix_chain_multiplication <- function(matrices, split = NULL) {
  n <- length(matrices)

  if(is.null(split)){
    split <- get_split_mat(matrices)
  }

  # Perform matrix multiplication in optimal order using split matrix
  multiply_matrices <- function(i, j) {
    if (i == j) {
      return(matrices[[i]])
    } else {
      k <- split[i, j]
      left <- multiply_matrices(i, k)
      right <- multiply_matrices(k + 1, j)
      return(mat_mult(left, right))
    }
  }

  result_matrix <- multiply_matrices(1, n)

  return(result_matrix)
}
