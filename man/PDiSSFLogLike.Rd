% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PDiSSFLogLike.R
\name{PDiSSFLogLike}
\alias{PDiSSFLogLike}
\title{PDiSSFLogLike}
\usage{
PDiSSFLogLike(beta, X1, X2, locations, k1, k2, maximumGap)
}
\arguments{
\item{beta}{Coefficients to be estimated}

\item{X1}{Data from the covariates for the iSSF component}

\item{X2}{Data from the covariates for the probability of detection component}

\item{locations}{data frame of GPS fixes}

\item{k1}{Number of covariates for the iSSF}

\item{k2}{Number of covariates for probability of detection. Equal to 0
if no probability of detection is being estimated.}

\item{maximumGap}{Maximum allowable number of consecutive missing fixes. Default is 3
(3 consecutive missing locations). Larger max lags will increase computing time
and may result in non convergence.}
}
\description{
Calculates the loglikelihood value for a set of parameters.
Not recommended for isolated use.
}
