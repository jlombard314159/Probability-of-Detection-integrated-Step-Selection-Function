\name{habitat}
\alias{habitat}
\docType{data}
\title{
Habitat covariates
}
\description{
A data frame of habitat covariates for all available cells in a study area.
}
\usage{data("habitat")}
\format{
  A data frame with 185 observations and 5 variables:
  \describe{
    \item{\code{utmX}}{a numeric vector of X coordinate}
    \item{\code{utmY}}{a numeric vector of Y coordinate}
    \item{\code{unitID}}{a numeric vector of cell IDs}
    \item{\code{prctSage}}{a numeric vector of Percent Big-sage}
    \item{\code{elevation}}{a numeric vector of elevation (km)}
  }
}
\source{
Nielson, R. M., B. F. J. Manly, L. L. McDonald, H. Sawyer, and T. L. McDonald. 2009. Estimating habitat selection when GPS fix success is less than 1. Ecology 90:2956-2962.
}
\examples{
data("habitat")
str(habitat)
}
\keyword{datasets}
