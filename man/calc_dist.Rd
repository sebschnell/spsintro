% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calc_dist}
\alias{calc_dist}
\title{Calculate the distance between two observations based on a set of auxiliary
variables}
\usage{
calc_dist(x_i, x_j, M)
}
\arguments{
\item{x_i}{Vector of auxiliary variables for observation \code{i}}

\item{x_j}{Vector of auxiliary variables for observation \code{j}}

\item{M}{A square matrix. When \code{M} is the identity matrix the
Euclidean distance results.}
}
\description{
Calculate the distance between two observations based on a set of auxiliary
variables
}
