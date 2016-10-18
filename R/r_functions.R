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
