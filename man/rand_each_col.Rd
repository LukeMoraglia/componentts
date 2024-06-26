% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/randomize_mat.R
\name{rand_each_col}
\alias{rand_each_col}
\title{Randomize (permute) each column of a matrix}
\usage{
rand_each_col(X)
}
\arguments{
\item{X}{the matrix to permute}
}
\value{
a matrix with the values in each column of \code{X} reordered
}
\description{
This function randomizes or permutes a matrix by reordering the values in each
column of the matrix. This results in a matrix where the relationships (correlations)
within the matrix are broken, and when compared to another matrix with different
variables and the same participants, the relationships between the two matrices
are also now broken.
}
\details{
This function is the option used when selecting "rand_each_col" in the stopping rule
functions such as \code{\link[=SRP]{SRP()}} or \code{\link[=DRP]{DRP()}}. However, all the stopping rule functions
in this package use a faster version written with \code{Rcpp}.
}
