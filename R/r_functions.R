#' A function to calculate fibonacci numbers
#'
#' @param n The length of the sequence
#'
#' @return A vector containing \code{n} fibonacci numbers
#' @export
#'
#' @examples
#' n <- 20;
#' gen_fib_r(n);
gen_fib_r <- function(n) {
  res <- rep(1, n);
  if (n <= 0) {
    return(0);
  } else if (n <= 2) {
    return(res);
  } else {
    a <- 1;
    b <- 1;
    for (i in seq(3, n, 1)) {
      c <- a + b;
      res[i] <- c;
      if (n <= 3) break();
      a <- b;
      b <- c;
    }
  }
  return(res);
}

#' @export
add_x <- function(y) {
  return(y + x);
}


#' The factorial of an integer number
#'
#' Meant to use only for teaching purposes not for replacing the built-in
#' function \code{factorial()}.
#' @param n An integer number
#'
#' @details A base R implementation using a \code{for} loop.
#'
#' @return The factorial of \code{n}
#' @export
#'
#' @examples
#' fact_r(10);
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


#' The factorial of an integer number
#'
#' Meant to use only for teaching purposes not for replacing the built-in
#' function \code{factorial()}.
#' @param n An integer number
#'
#' @details A base R implementation using recursive programming.
#'
#' @return The factorial of \code{n}
#' @export
#'
#' @examples
#' fact_rec_r(10);
fact_rec_r <- function(n) {
  if (n < 0) {
    stop("No factorial for negative numbers");
  } else if (n == 0) {
    return(1);
  } else {
    return(n*Recall(n - 1));
  }
}
