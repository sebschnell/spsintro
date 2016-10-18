#include "utilities.h"

//' Return the area of an irregular polygon
//' @param x A vector of x-coordinates
//' @param y A vector of y-coordinates (positive \code{y}-vaues mean upwards)
//' @return The area of the polygon
//' @export
// [[Rcpp::export]]
double calc_poly_area(const std::vector<double> &x,const std::vector<double> &y) {
  double a = 0.0;
  size_t n = x.size();
  auto j = n - 1; // the last element is the previous to the first
  for (size_t i = 0; i != n; ++i) {
    a += (x[j] + x[i]) * (y[j] - y[i]);
    j = i; //j is previous element to i
  }
  return abs(a)/2;
}
