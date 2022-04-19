# Probability-of-Detection-integrated-Step-Selection-Function (PDiSSF)

This R package comes with no guarantees, and technical support is limited. See Supporting Information from Vales et al. (2022; link here).  

We appreciate thoughtful suggestions for future improvements. Future development of this package may include using raster files for covariate layers and simplifying data preparation for the user.

# Installing the Package

The recommended approach for installing PDiSSF is to use download the .zip file and install to R locally. 
The recommended specifications for installing the package are: R  version blah blah.



## This package estimates the probability of detection integrated step selection function (PDiSSF) for animal movement data collected using global positioning system (GPS) units. 

GPS collars or backpacks attached to animals can fail to connect to satellites for various reasons, including habitat characteristics. Suppose GPS fix attempts fail when the animal is under thick canopy cover, on a steep slope, or somewhere that visible sky is restricted. In that case, the resulting location data could be biased away from these types of habitats. Analysis of habitat (resource) selection by animals using partial location data can result in biased conclusions unless fix success is accounted for in the study. This includes the estimation of step selection functions (Avgar et al. 2016).

This package estimates the probability of detection integrated step selection function (PDiSSF). The PDiSSF was developed to simultaneously model animal movement and habitat selection and the probability of detection given habitat conditions. The PDiSSF is described in Vales et al. (2022) but was initially presented as the "RSF for GPS fix success" in Nielson et al. (2009).

### Authors 
John Lombardi (jlombard314@gmail.com) and Ryan Nielson (Eagle Environmental, Inc.)

Maintainer: jlombard314@gmail.com

Other contributors: David Vales (Muckleshoot Indian Tribe), Michael Middleton (Muckleshoot Indian Tribe), and 
  Trent McDonald (McDonald Data Sciences, LLC)


#### References

Avgar, T., J. R. Potts, M. A. Lewis, and M. S. Boyce. 2016. Integrated step selection analysis: bridging the gap between resource selection and animal movement. Methods in Ecology and Evolution 7:619–630.

Nielson, R. M., B. F. J. Manly, L. L. McDonald, H. Sawyer, and T. L. McDonald. 2009. Estimating habitat selection when GPS fix success is less than 1. Ecology 90:2956-2962.

Vales, D, R. M. Nielson, and M. P. Middleton. 2022. Black-tailed Deer Seasonal Habitat Selection: Accounting for Missing GPS Fixes. Journal Vol:page-page.

