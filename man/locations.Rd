\name{locations}
\alias{locations}
\docType{data}
\title{
Animal locations (cellIDs)
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("locations")}
\format{
  A data frame with 444 observations on the following 2 variables.
  \describe{
    \item{\code{unitID}}{a numeric vector of cell IDs}
    \item{\code{fix}}{a numeric vector of fix attempt IDs}
  }
}
\details{
cell ID values should represent which cell the animal was in at every fix attempt, so NA's exist when the animal
wasn't located. The 'fix' column is not necessary for analysis, but hopefully reminds the 
user that the location attempts should be in sequential order and that unsuccessful 
fix attempts need to be represented in the data.
}
\source{
Nielson, R. M., B. F. J. Manly, L. L. McDonald, H. Sawyer, and T. L. McDonald. 2009. Estimating habitat selection when GPS fix success is less than 1. Ecology 90:2956-2962.
}
\examples{
data("locations")
str(locations)
}
\keyword{datasets}
