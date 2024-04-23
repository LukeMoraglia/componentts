#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

//' Fast matrix multiplication using C++
//'
//' This function offers a faster matrix multiplication than R provides natively.
//' It takes two matrices, `A` and `B` and returns the matrix multiplication of
//' `A %*% B`.
//'
//' @param A left matrix in the multiplication
//' @param B right matrix in the multiplication
//'
//' @returns a matrix of size `nrow(A)` by `ncol(B)`
//' @export
//'
// [[Rcpp::export]]
SEXP mat_mult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;

  return Rcpp::wrap(C);
}
