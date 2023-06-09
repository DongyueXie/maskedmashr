% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maskedmashFDR.R
\name{maskedmashFDR_condition}
\alias{maskedmashFDR_condition}
\title{masked mash for FDR control for each condition}
\usage{
maskedmashFDR_condition(obj, alpha = 0.05)
}
\arguments{
\item{obj}{fitted maskedmash object, or a list of P, result$lfsr, p.thresh}

\item{alpha}{target FDR level}
}
\value{
a list of rej.set indx and fdp_hat
}
\description{
order the tests from most to least significant based on lfsr from masked.mash, and remove least significant ones until FDP<alpha
}
