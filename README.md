# Probability-of-Detection-integrated-Step-Selection-Function (PDiSSF)

This R package comes with no guarantees, and technical support is limited. See Supporting Information from Vales et al. (2022; link here).  

We appreciate thoughtful suggestions for future improvements. Future development of this package may include using raster files for covariate layers and simplifying data preparation for the user.

# Installing the Package

The recommended approach for installing PDiSSF is to use download the .zip file and install to R locally. 
The recommended specifications for installing the package are: R >= 4.0.0. 



## This package estimates the probability of detection integrated step selection function (PDiSSF) for animal movement data collected using global positioning system (GPS) units. 

GPS collars or backpacks attached to animals can fail to connect to satellites for various reasons, including habitat characteristics. Suppose GPS fix attempts fail when the animal is under thick canopy cover, on a steep slope, or somewhere that visible sky is restricted. In that case, the resulting location data could be biased away from these types of habitats. Analysis of habitat (resource) selection by animals using partial location data can result in biased conclusions unless fix success is accounted for in the study. This includes the estimation of step selection functions (Avgar et al. 2016).

This package estimates the probability of detection integrated step selection function (PDiSSF). The PDiSSF was developed to simultaneously model animal movement and habitat selection and the probability of detection given habitat conditions. The PDiSSF is described in Vales et al. (2022) but was initially presented as the "RSF for GPS fix success" in Nielson et al. (2009).

## Example data

The package provides some example data for running, along with examples in the specific help files.

Additional data are provided in the 'exampleHabitat' and 'exampleLocations' .RData files in the ./data folder.


#### exampleHabitat.RData
A data.frame named habitat containing covariate values for all 30 m cells within the deer’s home range. Data from the analysis by Vales et al. (2022). See the manuscript for details. 
This data.frame contains the following information:
pointID: unique ID for each 30 m cell in the deer’s home range
cc: Percent canopy cover for a cell
cc75: cancopy cover shifted (centered) by subtracting 75 meters. 
ab: accepted biomass (amount of neutral and selected forage)
dcfe: distance (m) to cover-forage edge
drds: distance (m) to nearest road
ealbers: easting coordinate for the cell in UTM
nalbers: northing coordinate for the cell in UTM

#### exampleLocations.RData
A data.frame named locations of GPS fixes from one deer used in the analysis by Vales et al. (2022). See the manuscript for details. This data.frame contains the following information:
lineno: Ignored. This is simply to recognize that all fix attempts (missing and successful) need to be represented in the data. The numbers should be consecutive and have no gaps.
pointID: cell ID for 30 m cells where the deer was located at each fix attempt. Missing locations should be in the data and represented by NA. All pointID values should be represented in the habitat data.frame. 

Example model run based on 'exampleHabitat' and 'exampleLocation' data sets. These
are provided in the ./data folder.

Timing is based on:
Lenovo Yoga
1.8 Ghz 
Intel i7
16gb RAM

``` library(PDiSSF)
#takes about 2 minutes
pdissf(habitatDF = habitat,
       CellID = locations$pointID,
       habitatCellID = 'pointID',
       iSSFCovars = c("ab", "dcfe", "drds", "distance"),
       probDetCovars = NULL,
       distColumns  = c("ealbers","nalbers"),
       maximumGap = 3
)

# fit PDiSSF
# takes ~50 minutes
pdissf(habitatDF = habitat,
       CellID = locations$pointID,
       habitatCellID = 'pointID',
       iSSFCovars = c("ab", "dcfe", "drds", "distance"),
       probDetCovars = c("cc75"),
       distColumns  = c("ealbers","nalbers"),
       maximumGap = 3
)

```

### Authors 
John Lombardi (jlombard314@gmail.com) and Ryan Nielson (Eagle Environmental, Inc.; ryan@eagleenvironmental.net)

Maintainer: jlombard314@gmail.com

Other contributors: David Vales (Muckleshoot Indian Tribe), Michael Middleton (Muckleshoot Indian Tribe), and 
  Trent McDonald (McDonald Data Sciences, LLC)


#### References

Avgar, T., J. R. Potts, M. A. Lewis, and M. S. Boyce. 2016. Integrated step selection analysis: bridging the gap between resource selection and animal movement. Methods in Ecology and Evolution 7:619–630.

Nielson, R. M., B. F. J. Manly, L. L. McDonald, H. Sawyer, and T. L. McDonald. 2009. Estimating habitat selection when GPS fix success is less than 1. Ecology 90:2956-2962.

Vales, D, R. M. Nielson, and M. P. Middleton. 2022. Black-tailed Deer Seasonal Habitat Selection: Accounting for Missing GPS Fixes. Journal Vol:page-page.

