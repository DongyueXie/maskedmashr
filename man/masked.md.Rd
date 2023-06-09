% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maskedmd.R
\name{masked.md}
\alias{masked.md}
\title{masked deconvolution: estimate prior covariance matrices using multivariate deconvolution and masked data}
\usage{
masked.md(
  Z,
  Z.comb,
  strong = NULL,
  thresh = NULL,
  Pi = NULL,
  U.canon = NULL,
  U.data = NULL,
  npc = 5,
  usepointmass = FALSE,
  algorithm.version = "Rcpp",
  adjust = "lb",
  nu = NULL,
  pi_thresh = 1e-08,
  max_iter = 1000,
  tol = 1e-05,
  verbose = TRUE,
  printevery = 50,
  adj.const = 0.001
)
}
\arguments{
\item{strong}{index of samples that have "strong" effects.}

\item{thresh}{NULL or a number, |z-score| larger than thresh or |z|<=s(t) will be masked. when thresh<=Phi^{-1}(0.75), all masked.}

\item{Pi}{prior weights; either fixed or use as init value}

\item{U.canon}{a list of canonical prior cov matrices; fixed.}

\item{U.data}{a list of data-driven prior cov matrices; If null, will use pca on masked z score for initialization.}

\item{algorithm.version}{either Rcpp or R for calculating likelihood}

\item{adjust}{adjust diagonal of estimated U}

\item{nu}{prior df of U. If NULL, nu = R+1, for adjust = 'prior'}

\item{data}{object from `mash_set_data` function}
}
\description{
This function estimates prior weights and variances using masked z scores. If strong
is provided, then U.canon should be set to NULL and usepointmass set to FALSE;
}
