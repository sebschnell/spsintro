library('Rcpp');
setwd("C:/Users/Schnell/Documents/Courses/introduction_programming_simulation");


#'Calculate the factorial of an integer number
#'
#'Meant to use only for teaching purposes not for replacing the built-in
#'function \code{factorial()}.
#'
#'@param n An integer number
#'
#'@details Internally a \code{for} loop is used.
#'
#'@return The factorial of \code{n}
#'@export
#'
#' @examples
#'    fact_r(10);
fact_r <- function(n) {
  f <- 1;
  if (n < 0) {
    stop("No factorial for negative numbers");
  } else if (n == 0) {
    return(1);
  } else {
    for (i in seq.int(n)) {
      f = f*i;
    }
    return(f);
  }
}

#'Calculate the factorial of an integer number
#'
#'Meant to use only for teaching purposes not for replacing the built-in
#'function \code{factorial()}.
#'
#'@param n An integer number
#'
#'@details This function works using recursion (the function calls itself).
#'
#'@return The factorial of \code{n}
#'@export
#'
#' @examples
#'    fact_rec_r(10);
fact_rec_r <- function(n) {
  if (n < 0) {
    stop("No factorial for negative numbers");
  } else if (n == 0) {
    return(1);
  } else {
    return(n*Recall(n - 1));
  }
}
