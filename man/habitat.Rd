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
    \item{\code{utmX}}{X coordinate in UTM}
    \item{\code{utmY}}{Y coordinate in UTM}
    \item{\code{unitID}}{Cell IDs. Must use the same column name in the "locations" data frame.}
    \item{\code{prctSage}}{Percent Big-sage} 
    \item{\code{elevation}}{Elevation (km)}
  }
  NOTE: missing values are not allowed.
}
\source{
Nielson, R. M., B. F. J. Manly, L. L. McDonald, H. Sawyer, and T. L. McDonald. 2009. Estimating habitat selection when GPS fix success is less than 1. Ecology 90:2956-2962.
}
\examples{
data("habitat")
str(habitat)
}
\keyword{datasets}
