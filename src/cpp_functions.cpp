#include <RcppEigen.h>
// Use quotation marks to search for header files that are located in the same
// directory as the file containing the directive
#include "Tree.h"
#include "Polygon.h"
#include "utilities.h"
#include <vector> // Brackets for implementation dependent search paths
#include <algorithm>
#include <utility>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]



//' A function to calculate fibonacci numbers (implemented in C++)
//'
//' @param n The length of the sequence
//' @return A vector containing \code{n} fibonacci numbers
//' @examples
//' n <- 20;
//' gen_fib_cpp(n);
//' @export
// [[Rcpp::export]]
std::vector<double> gen_fib_cpp(const int n) {
  if (n <= 0) {
    std::vector<double> res(1, 0);
    return res;
  } else if (n <= 2) {
    std::vector<double> res(n, 1);
    return res;
  } else {
    std::vector<double> res(n, 1);
    double a = 1, b = 1;
    double c;
    for (int i = 2; i < n; ++i) {
      c = a + b;
      res[i] = c;
      if (n <= 3) break;
      a = b;
      b = c;
    }
    return res;
  }
}

//' Find the maximum from a list of values
//'
//' @param x A vector holding a number of values
//' @return The largest number in \code{x}
//' @export
// [[Rcpp::export]]
double find_max(const std::vector<double> &x) {
  double temp = 0.0;
  for (unsigned i = 0; i != x.size(); ++i) {
    if (x[i] > temp) temp = x[i];
  }
  return(temp);
}

//' Calculate the distance between two observations based on a set of auxiliary
//' variables
//'
//' @param x_i Vector of auxiliary variables for observation \code{i}
//' @param x_j Vector of auxiliary variables for observation \code{j}
//' @param M A square matrix. When \code{M} is the identity matrix the
//' Euclidean distance results.
//' @export
// [[Rcpp::export]]
double calc_dist(const Eigen::VectorXd & x_i,
                 const Eigen::VectorXd & x_j,
                 const Eigen::MatrixXd & M) {
  return sqrt((x_i - x_j).transpose()*M*(x_i - x_j));
}


//' Search nearest neighbours
//'
//' @param X Matrix with auxiliary variables of the reference set. Rows are the
//' single observations and variables are stored in columns.
//' @param M A square matrix defining the distance measure. When \code{M} is
//' the identity matrix, the Euclidean distance results.
//' @param k The number of neighbours
//' @export
// [[Rcpp::export]]
Rcpp::List find_knn(const Eigen::MatrixXd & Y,
                    const Eigen::MatrixXd & X,
                    const Eigen::MatrixXd & M,
                    const unsigned & k) {
  Eigen::MatrixXi NN(Y.rows(), k);
  Eigen::MatrixXd D(Y.rows(), k);

  std::vector<std::pair<double, unsigned> > d;
  d.reserve(Y.rows());
  unsigned idx;
  for (unsigned i = 0; i < Y.rows(); ++i) {
    d.clear();
    for (unsigned j = 0; j < X.rows(); ++j) {
      ++idx;
      if (i != j) {
        d.push_back({calc_dist(Y.row(i), X.row(j), M), (j + 1)});
      }
    }
    std::partial_sort(d.begin(), d.begin() + k, d.end());
    for (unsigned k_i = 0; k_i < k; ++k_i) {
      NN(i, k_i) = d[k_i].second;
      D(i, k_i) = d[k_i].first;

    }
  }
  int n_var = Y.cols();
  return Rcpp::List::create (Rcpp::Named("nn") = NN,
                             Rcpp::Named("dist") = D,
                             Rcpp::Named("k") = k,
                             Rcpp::Named("n_var") = n_var);
}

//' Print C++ data type ranges to the R-console
//' @return No return values, simply printing to the R-console
//' @export
// [[Rcpp::export]]
void show_type_ranges() {
  Rcpp::Rcout << "bool " << std::numeric_limits<bool>::min() << "-" << std::numeric_limits<bool>::max() << ", " << sizeof(bool) << " byte" << std::endl;
  int minc = std::numeric_limits<char>::min(); // Convert to integer to see numbers not character values
  int maxc = std::numeric_limits<char>::max();
  Rcpp::Rcout << "char " << minc << "-" << maxc << ", " << sizeof(char) << " byte" << std::endl;
  Rcpp::Rcout << "short " << std::numeric_limits<short>::min() << "-" << std::numeric_limits<short>::max() << ", " << sizeof(short) << " bytes" << std::endl;
  Rcpp::Rcout << "int " << std::numeric_limits<int>::min() << "-" << std::numeric_limits<int>::max() << ", " << sizeof(int) << " bytes" << std::endl;
  Rcpp::Rcout << "unsigned int " << std::numeric_limits<unsigned int>::min() << "-" << std::numeric_limits<unsigned int>::max() << ", " << sizeof(unsigned int) << " bytes" << std::endl;
  Rcpp::Rcout << "long " << std::numeric_limits<long>::min() << "-" << std::numeric_limits<long>::max() << ", " << sizeof(long) << " bytes" << std::endl;
  Rcpp::Rcout << "long long " << std::numeric_limits<long long>::min() << "-" << std::numeric_limits<long long>::max() << ", " << sizeof(long long) << " bytes" << std::endl;
  Rcpp::Rcout << "float " << std::numeric_limits<float>::min() << "-" << std::numeric_limits<float>::max() << ", " << sizeof(float) << " bytes" << std::endl;
  Rcpp::Rcout << "double " << std::numeric_limits<double>::min() << "-" << std::numeric_limits<double>::max() << ", " << sizeof(double) << " bytes" << std::endl;
  Rcpp::Rcout << "long double " << std::numeric_limits<long double>::min() << "-" << std::numeric_limits<long double>::max() << ", " << sizeof(long double) << " bytes" << std::endl;
}

//' Print C++ data type conversions
//' @return No return values, simply printing to the R-console
//' @export
// [[Rcpp::export]]
void type_conversions() {
  bool b = 42;
  Rcpp::Rcout << b << std::endl;
  int i = b;
  Rcpp::Rcout << i << std::endl;
  i = 3.14;
  Rcpp::Rcout << i << std::endl;
  double pi = i;
  Rcpp::Rcout << pi << std::endl;
  unsigned char c = -1;
  int cv = c;
  Rcpp::Rcout << cv << std::endl;
  signed char c2 = 256;
  Rcpp::Rcout << c2 << std::endl;
}

//' Show effect of mixing signed and unsigned types
//' @return No return values, simply printing to the R-console
//' @export
// [[Rcpp::export]]
void mix_signed_unsigned() {
  int a = -1, b = 1;
  unsigned bu = 1;
  Rcpp::Rcout << "int*int " << a*b << std::endl;
  Rcpp::Rcout << "int*unsigned " << a*bu << std::endl;
}

//' Illustrate scoping rules
//' @return No return values, simply printing to the R-console
//' @export
// [[Rcpp::export]]
void illustrate_scope() {
  int i = 100, sum = 0;
  for (int i = 0; i != 10; ++i) {
    sum += i;
  }
  Rcpp::Rcout << i << " " << sum << std::endl;
}

//' Function to illustrate references
//' @return No return values, simply printing to the R-console
//' @export
// [[Rcpp::export]]
void show_references() {
  int i = 256;
  Rcpp::Rcout << "i: " << i << std::endl;
  int &ref_val = i;
  Rcpp::Rcout << "ref_val: " << ref_val << std::endl;
  ref_val = 2;
  Rcpp::Rcout << "ref_val = 2: " << i << std::endl;
  int &ref_val2 = ref_val;
  Rcpp::Rcout << "ref_val2: " << ref_val2 << std::endl;
  int ii = ref_val;
  Rcpp::Rcout << "ii: " << ii << std::endl;
}

//' Function to illustrate pointers
//' @return No return values, simply printing to the R-console
//' @export
// [[Rcpp::export]]
void show_pointers() {
  double dval = 3.14;
  double *pd = &dval;
  Rcpp::Rcout << "dval: " << dval << std::endl;
  Rcpp::Rcout << "pd: " << pd << std::endl;
  Rcpp::Rcout << "*pd: " << *pd << std::endl;
  *pd = 0.0;
  Rcpp::Rcout << "*pd = 0.0: " << dval << std::endl;
}

//' Illustrates definition and usage of a simple class with only data members
//' @param dbh The trees diameter at breast height (double)
//' @param height The trees height (double)
//' @param no_stems The number of stems the tree has (int)
//' @param spec_name The species of the tree (string)
//' @return No return values, simply printing to the R-console
//' @export
// [[Rcpp::export]]
void create_tree_struct(double dbh, double height, int no_stems, std::string spec_name) {
  Tree tree1;
  tree1.dbh = dbh;
  tree1.height = height;
  tree1.no_stems = no_stems;
  tree1.spec_name = spec_name;
  Rcpp::Rcout << tree1.spec_name << std::endl;
  Rcpp::Rcout << &tree1 << std::endl;
  Rcpp::Rcout << sizeof(tree1) << std::endl;
}

//' Convert all characters in a string to upper case
//' @param str A string
//' @return The same string as \code{str} but with upper case letters
//' @export
// [[Rcpp::export]]
std::string to_upper_range_for(std::string str) {
  for (auto &c : str)
    c = toupper(c);
  return str;
}

//' Convert all characters in a string to upper case
//' @param str A string
//' @return The same string as \code{str} but with upper case letters
//' @export
// [[Rcpp::export]]
std::string to_upper_subscript(std::string str) {
  for (decltype(str.size()) index = 0; index != str.size(); ++index)
    str[index] = toupper(str[index]);
  return str;
}

//' Add elements to a vector
//' @param n An integer giving the final size of the vector
//' @return A vector of size \code{n}
//' @export
// [[Rcpp::export]]
std::vector<int> add_elem(int n) {
  std::vector<int> v;
  for (int i = 0; i != n; ++i)
    v.push_back(i);
  return v;
}

//' Sums all elements in a vector
//' @param v A vector of doubles
//' @return A double containing the sum of all elements of \code{v}
//' @export
// [[Rcpp::export]]
double sum_elem(std::vector<double> v) {
  double sum = 0.0;
  for (auto &i : v)
    sum += i;
  return sum;
}

//' Square vector elements
//' @param v A vector of doubles
//' @return A double containing the sum of all elements of \code{v}
//' @export
// [[Rcpp::export]]
std::vector<double> squ_elem(std::vector<double> v) {
  for (decltype(v.size()) i = 0; i != v.size(); ++i)
    v[i] *= v[i];
  return v;
}

//' Square vector elements
//' @param v A vector of doubles
//' @return A double containing the sum of all elements of \code{v}
//' @export
// [[Rcpp::export]]
std::vector<double> squ_elem_it(std::vector<double> v) {
  for (auto it = v.begin(); it != v.end(); ++it)
    *it *= *it;
  return v;
}

//' Array element access example
//' @return Nothing is returned, only printing to the console
//' @export
// [[Rcpp::export]]
void arr_access() {
  char c[] = "Hello World!";
  for (auto i : c)
    Rcpp::Rcout << i << "";
  Rcpp::Rcout << std::endl;

  auto sz = sizeof(c)/sizeof(*c);
  for (size_t i = 0; i != sz; ++i)
    Rcpp::Rcout << c[i] << "";
  Rcpp::Rcout << std::endl;
}

//' Return remainder
//' @return The remainder of an integer division
//' @export
// [[Rcpp::export]]
int ret_rem(int a, int b) {
  Rcpp::Rcout << "Integer division quotient: " << a/b << std::endl;
  Rcpp::Rcout << "Integer division remainder: " << a%b << std::endl;
  return a%b;
}

//' Conditional operator
//' @param x An int
//' @return Nothing to return, just printing to the console
//' @export
// [[Rcpp::export]]
void dem_cond(int x) {
  std::string win_loose =  x < 50 ? "win" : "loose";
  Rcpp::Rcout << win_loose << std::endl;
}

//' An if-else example
//' @param dbh Diameter at breast height in cm
//' @return A diameter class
//' @export
// [[Rcpp::export]]
int dem_if_else(double dbh) {
  int dbh_class = 0;
  if (dbh > 50.0)
    dbh_class = 3;
  else {
    if (dbh <= 20)
      dbh_class = 1;
    else
      dbh_class = 2;
  }
  return dbh_class;
}

//' Find first negative element in a vector and return its position
//' @param x A vector of double
//' @return Index of the first negative value
//' @export
// [[Rcpp::export]]
int find_first_neg(std::vector<double> x) {
  int idx = 0;
  auto it = x.begin();
  while (it != x.end() && *it >= 0) {
    ++it;
    ++idx;
  }
  if (it == x.end())
    idx = 0;
  return idx + 1; // In R indexes begin with 1
}

//' Find first negative element in a vector and return its position
//' @param x A vector of double
//' @return Index of the first negative value
//' @export
// [[Rcpp::export]]
int find_first_neg_for(std::vector<double> x) {
  int idx = 0;
  auto it = x.begin();
  for (; it != x.end() && *it >= 0; ++it) {
    ++idx;
  }
  if (it == x.end())
    idx = 0;
  return idx + 1; // In R indexes begin with 1
}

//' Find first negative element in a vector and return its position
//' @param x A vector of double
//' @return Index of the first negative value
//' @export
// [[Rcpp::export]]
int find_first_neg_rf(std::vector<double> x) {
  int idx = 0;
  for (auto &r : x) {
    if (r < 0)
      break;
    ++idx;
  }
  return idx + 1; // In R indexes begin with 1
}

//' Find first negative element in a vector and return its position
//' @param x A vector of double
//' @return Index of the first negative value
//' @export
// [[Rcpp::export]]
unsigned find_first_neg_dw(std::vector<double> x) {
  unsigned idx = 0;
  do {
    ++idx;
  } while ((idx - 1) != x.size() && x[idx - 1] >= 0);

  if ((idx - 1) == x.size())
    idx = 0;
  return idx;
}

//' Demonstrate rounding in C++
//' R's built-in \code{round()} function shows a
//' different behaviour than the \code{round()} function of C++
//' @param x A double value to be rounded
//' @return The rounded integer represented as double
//' @export
// [[Rcpp::export]]
double round_cpp(double x) {
  return round(x);
}

//' Test Polygon class
//' @return Return nothing
//' @export
// [[Rcpp::export]]
void use_polygon_class() {
  std::vector<double> x{0, 6, 6, 0};
  std::vector<double> y{0, 0, 4, 4};
  Rectangle rect(x, y, 6, 4);
  Rcpp::Rcout << "height: " << rect.get_h() << std::endl;
  Rcpp::Rcout << "width " << rect.get_w() << std::endl;
  Rcpp::Rcout << "area: " << rect.area() << std::endl;
}

//' Demonstrates inheritance and dynamic binding
//' @return Return nothing
//' @export
// [[Rcpp::export]]
void inheritance() {
  std::vector<double> x_poly{0, 4, 4, 2, 0};
  std::vector<double> y_poly{0, 0, 4, 6, 4};
  std::vector<double> x_rect{0, 4, 4, 0};
  std::vector<double> y_rect{0, 0, 4, 4};
  std::vector<double> x_tri{0, 4, 2};
  std::vector<double> y_tri{4, 4, 6};

  Polygon poly(x_poly, y_poly);
  Rectangle rect(x_rect, y_rect, 4, 4);
  Triangle tri(x_tri, y_tri, 4, 2);

  Rcpp::Rcout << "Area poly: " << poly.area() << std::endl;
  Rcpp::Rcout << "Area rect: " << rect.area() << std::endl;
  Rcpp::Rcout << "Area tri: " << tri.area() << std::endl;
}

//' Naive implementation of bubble sort algorithm
//' @param x A numeric vector which values are sorted in increasing order
//' @return A numeric vector with values sorted in increasing order
// [[Rcpp::export]]
std::vector<double> bubble_sort (std::vector<double> &x) {
  double temp = 0.0;
  for (unsigned i = 0; i != x.size(); ++i) {
    for (unsigned j = 0; j != x.size() - 1; ++j) {
      if (x[j] > x[j + 1]) {
        temp = x[j];
        x[j] = x[j + 1];
        x[j + 1] = temp;
      }
    }
  }
  return(x);
}



