X <- matrix(1:30, 10, 3)
X

set.seed(42)
componentts:::rand_each_col_cpp(X)

set.seed(42)
componentts:::rand_rows_cpp(X)

