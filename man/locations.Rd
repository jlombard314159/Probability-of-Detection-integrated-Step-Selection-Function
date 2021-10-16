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
    \item{\code{fix}}{a vector of integers representing sequential fix attempts. Locations are sorted accordingly.
          fix = 1 for first attempt (always successful)
          fix = 2 for second attempt (successful or not)
          fix = n for the nth (final) fix attempt (always successful)}
  }
}
\details{
cell ID values should represent which cell the animal was in at every fix attempt, so NA's exist when the animal
wasn't located. The 'fix' column is not necessary for analysis, but hopefully reminds the 
user that the location attempts should be in sequential order and that unsuccessful 
fix attempts need to be represented in the data.

This example fits the PDiSSF described in Nielson et al. (2009) to data from that same paper. 

The computing demand can be high for larger data sets, depending fix success, the number of available cells, the number of fix attempts,and the number of covariates in the model. Fitting the PDiSSF could take seconds or days. The size of the study area (number of cells) has the most influence on memory requirements, processing time, and convergence. Parallel processing (try library parallel) can help when analyzing multiple animals, multitheading with Mircrosf R Open can speed up procesing. Parallel processing will not speed up model estimation, but merely allow the user to fit more than one model at a time. See the Supporting Information in Vales et al. (2022; __github__ for more information). 
}
\source{
Nielson, R. M., B. F. J. Manly, L. L. McDonald, H. Sawyer, and T. L. McDonald. 2009. Estimating habitat selection when GPS fix success is less than 1. Ecology 90:2956-2962.
}
\examples{
data("locations")
str(locations)
}
\keyword{datasets}
