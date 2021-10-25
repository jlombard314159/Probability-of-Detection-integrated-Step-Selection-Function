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
    \item{\code{unitID}}{Cell IDs. Must use the same name in the "habitat" file.}
    \item{\code{fixAttempt}}{Integers representing sequential fix attempts. Locations are sorted accordingly.
          fix = 1 for first attempt (always successful)
          fix = 2 for second attempt (successful or not)
          fix = n for the nth (final) fix attempt (always successful)}
  }
}
\details{
Cell ID values should represent which cell the animal was in at every fix attempt, so NA's exist when the animal
wasn't located. The 'fix' column is not necessary for analysis, because we have assumed that you have sorted the data by 
fix attempt. Thus, this is just a reminder to double-check. 

This data is from Nielson et al. (2009). 

}
\source{
Nielson, R. M., B. F. J. Manly, L. L. McDonald, H. Sawyer, and T. L. McDonald. 2009. Estimating habitat selection when GPS fix success is less than 1. Ecology 90:2956-2962.
}
\examples{
data("locations")
str(locations)
}
\keyword{datasets}
