library('Rcpp');
setwd("C:/Users/Schnell/Documents/Courses/introduction_programming_simulation");

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

fact_rec_r <- function(n) {
  if (n < 0) {
    stop("No factorial for negative numbers");
  } else if (n == 0) {
    return(1);
  } else {
    return(n*Recall(n - 1));
  }
}

sourceCpp(file = "R/factorial_example1.cpp");

x <- 10;
fact_r(x);
fact_rec_r(x);
fact_cpp(-1);
fact_rec_cpp(-1);

library(microbenchmark);
microbenchmark(fact_r(x),
               fact_rec_r(x),
               fact_cpp(x),
               fact_rec_cpp(x),
               factorial(x));

