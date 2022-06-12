\name{exampleHabitat}
\alias{exampleHabitat}
\docType{data}
\title{
Habitat covariates
}
\description{
A data.frame named habitat containing covariate values for all 
30 m cells within the deer’s home range. Data from the analysis by Vales et al. (2022). 
see the manuscript for details. 
}
\usage{data("exampleHabitat")}
\format{
  \describe{
    \item{\code{pointID}}{unique ID for each 30 m cell in the deer’s home range}
    \item{\code{cc}}{Percent canopy cover for a cell}
    \item{\code{cc75}}{cancopy cover shifted (centered) by subtracting 75 meters. }
    \item{\code{ab}}{accepted biomass (amount of neutral and selected forage)} 
    \item{\code{dcfe}}{distance (m) to cover-forage edge}
    \item{\code{drds}}{distance (m) to nearest road}
    \item{\code{ealbers}}{easting coordinate for the cell in UTM}
    \item{\code{nalbers}}{northing coordinate for the cell in UTM}
  }

}
\source{
Vales et al. (2022)
}
\examples{
data("exampleHabitat")
str(exampleHabitat)
}
\keyword{datasets}
