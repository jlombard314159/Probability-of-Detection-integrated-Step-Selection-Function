\name{PDiSSF-package}
\alias{PDiSSF-package}
\alias{PDiSSF}
\docType{package}
\title{
This package estimates the probability of detection integrated step selection function (PDiSSF) for animal movement data collected using global positioning system (GPS) units. 
}
\description{

GPS collars or backpacks attached to animals can fail to connect to satellites for various reasons, including habitat characteristics. Suppose GPS fix attempts fail when the animal is under thick canopy cover, on a steep slope, or somewhere that visible sky is restricted. In that case, the resulting location data could be biased away from these types of habitats. Analysis of habitat (resource) selection by animals using partial location data can result in biased conclusions unless fix success is accounted for in the study. This includes the estimation of step selection functions (Avgar et al. 2016).

This package estimates the probability of detection integrated step selection function (PDiSSF). The PDiSSF was developed to simultaneously model animal movement and habitat selection and the probability of detection given habitat conditions. The PDiSSF is described in Vales et al. (2022) but was initially presented as the "RSF for GPS fix success" in Nielson et al. (2009).

}

\author{

John Lombardi (jlombard314@gmail.com) and Ryan Nielson (Eagle Environmental, Inc.)

Maintainer: jlombard314@gmail.com

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
