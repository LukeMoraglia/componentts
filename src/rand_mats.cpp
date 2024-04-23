#include <Rcpp.h>
using namespace Rcpp;

// Thanks Dirk! https://gallery.rcpp.org/articles/stl-random-shuffle/

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

//' Randomize (permute) each column of a matrix
//'
//' The C++ implementation gives this function more speed than the `rand_each_col()`
//' function written in just in R.
//' This function randomizes or permutes a matrix by reordering the values in each
//' column of the matrix. This results in a matrix where the relationships (correlations)
//' within the matrix are broken, and when compared to another matrix with different
//' variables and the same participants, the relationships between the two matrices
//' are also now broken.
//'
//' This function is the option used when selecting "rand_each_col" in the stopping rule
//' functions such as [SRP()] or [DRP()].
//'
//' @param X the matrix to permute
//'
//' @returns a matrix with the values in each column of `X` reordered
//' @export
// [[Rcpp::export]]
NumericMatrix rand_each_col_cpp(NumericMatrix X) {
  int nrow = X.nrow();
  int ncol = X.ncol();
  NumericMatrix result(nrow, ncol);

  for (int j = 0; j < ncol; ++j) {
    IntegerVector indices = Rcpp::seq(0, nrow - 1);
    std::random_shuffle(indices.begin(), indices.end(), randWrapper);
    for (int i = 0; i < nrow; ++i) {
      result(i, j) = X(indices[i], j);
    }
  }

  return result;
}

//' Randomize (permute) the rows of a matrix
//'
//' The C++ implementation gives this function more speed than the `rand_rows()`
//' function written in just in R.
//' This function randomizes or permutes a matrix by reordering complete rows of
//' the matrix. This results in a matrix where the relationships (correlations)
//' within the matrix are preserved. But when compared to another matrix with different
//' variables and the same participants, the relationships between the two matrices
//' are now broken.
//'
//' This function is the option used when selecting "rand_rows" in the stopping rule
//' functions such as [SRP()] or [DRP()].
//'
//' @param X the matrix to permute
//'
//' @returns a matrix with the rows of `X` reordered
//' @export
// [[Rcpp::export]]
NumericMatrix rand_rows_cpp(NumericMatrix X){
  int nrow = X.nrow();
  int ncol = X.ncol();
  NumericMatrix result(nrow, ncol);

  IntegerVector indices = Rcpp::seq(0, nrow - 1);
  std::random_shuffle(indices.begin(), indices.end(), randWrapper);

  for (int j = 0; j < ncol; ++j) {
    for (int i = 0; i < nrow; ++i) {
      result(i, j) = X(indices[i], j);
    }
  }

  return result;

}
