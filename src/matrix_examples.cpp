#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]


//'Demonstrate the use of a numeric matrix object from the Rcpp package
//'
//'@param x A numeric matrix
//'
//'@return A new matrix with the same dimensions as \code{x} filled with zeroes
//'@export
// [[Rcpp::export]]
Rcpp::NumericMatrix use_matrix(const Rcpp::NumericMatrix &x) {
  Rcpp::NumericMatrix m1(x);
  Rcpp::NumericMatrix m2(5, 4);

  Rcpp::Rcout << m1(1, 2) << std::endl;

  return m2;
}

//'Demonstrate the use of a list object from the Rcpp package
//'
//'@param x A numeric matrix
//'@param y A numeric vector
//'
//'@return A named list containing both \code{x} and \code{y}.
//'@export
// [[Rcpp::export]]
Rcpp::List use_list(const Rcpp::NumericMatrix &x, const Rcpp::NumericVector &y) {
  return Rcpp::List::create(Rcpp::Named("m") = x,
                            Rcpp::Named("v") = y);
}
