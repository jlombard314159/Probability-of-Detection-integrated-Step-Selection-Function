\name{exampleLocations}
\alias{exampleLocations}
\docType{data}
\title{
Location information
}
\description{
A data.frame named locations of GPS fixes from one deer used in 
the analysis by Vales et al. (2022). See the manuscript for details. 
}
\usage{data("exampleLocations")}
\format{
  \describe{
    \item{\code{lineno}}{Ignored. This is simply to recognize that all fix attempts 
          (missing and successful) need to be represented in the data. }
    \item{\code{pointID}}{cell ID for 30 m cells where the deer was located at each fix attempt.
          Missing locations should be in the data and represented by NA. All pointID 
          values should be represented in the habitat data.frame. }
}
\source{
Vales et al. (2022)
}
\examples{
data("exampleLocations")
str(exampleLocations)
}
\keyword{datasets}
