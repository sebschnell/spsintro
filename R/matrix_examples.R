rm(list = ls());
setwd("C:\\Users\\Schnell\\Documents\\Courses\\introduction_programming_simulation");
library('Rcpp');

sourceCpp("R\\matrix_examples.cpp");
x <- matrix(rnorm(20, 5, 1), nrow = 5, ncol = 4);
y <- runif(100);

use_matrix(x);
use_list(x, y);
