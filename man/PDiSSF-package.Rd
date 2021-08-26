\name{PDiSSF-package}
\alias{PDiSSF-package}
\alias{PDiSSF}
\docType{package}
\title{
Fits the probability of detection integrated step selection function (PDiSSF) to animal movement data collected via global positioning system (GPS) units. 
}
\description{
GPS collars or backpacks attached to animals can fail to connect to satelites for a variety of reasons, including habitat characteristics. If GPS fix attempts fail when the animal is under thick canopy cover, on a steep slope, or somewhere that visible sky is restricted, the resulting location data could be biased away from these types of habitat. Analysis of habitat (resource) selection by animals using biased location data can result in biased conclusions, unless fix success is accounted for in the analysis. This includes estimation of step selection functions (Avgar et al. 2016).

This package estimates the probability of detection integrated sep selection function (PDiSSF), which was developed to simultaneously model animal movement and habitat selection, along with probability of detection given habitat conditions. The PDiSSF is described in Vales et al. (2021), but was originally presented as the "RSF for GPS fix success" in Nielson et al. (2009).
}
\author{
Ryan Nielson (Eeagle Environmental, Inc.) and John Lombardi (WEST, Inc.).

Maintainer: jlomardi@west-inc.com

Other contributors: David Vales (Muckleshoot Indian Tribe), Michael Middleton (Muckleshoot Indian Tribe), and 
  Trent McDonald (McDonald Data Sciences, LLC)
}
\references{
Avgar, T., J. R. Potts, M. A. Lewis, and M. S. Boyce. 2016. Integrated step selection analysis: bridging the gap between resource selection and animal movement. Methods in Ecology and Evolution 7:619–630.

Nielson, R. M., B. F. J. Manly, L. L. McDonald, H. Sawyer, and T. L. McDonald. 2009. Estimating habitat selection when GPS fix success is less than 1. Ecology 90:2956-2962.

Vales, D, R. M. Nielson, and M. P. Middleton. 2022. Black-tailed Deer Seasonal Habitat Selection: Accounting for Missing GPS Fixes. Journal Vol:page-page.


}
\keyword{ package }
\seealso{
\code{\link{pdissf}}
\code{\link{habitat}}
\code{\link{locations}}
}