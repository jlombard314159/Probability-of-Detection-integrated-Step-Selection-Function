\name{pdissf}
\alias{pdissf}
\title{
Fits the probability of detection integrated step selection function (PDiSSF)
}
\description{
This package fits the probability of detection integrated detection fucntion presented in Vales et al. (2001), and presented as the "RSF for GPS fix succes" in Nielson et al. (2009). As of now, this package relies on 
totally pre-processed data. See help(habitat) and help(locations) for more details and examples.
}
\usage{
pdissf(habitatDF = habitat, 
      CellID = locations$unitID, 
      iSSFCovars = c("distance", "prctSage", "elevation"), 
      probDetCovars = "prctSage", 
      distColumns = c("utmX","utmY"), 
      maxLagArg = 4)
      
This example fits the PDiSSF described in Nielson et al. (2009) to data from that same paper. 

Note that computing demand can be high, depending fix success, the number of available cells, the number of fix attempts,
and the number of covariates in the model. Fitting the PDiSSF could take seconds or days. Parallel processing can help with large batches, but...DAVE to finish this sentence. See Vales et al. (2021; __github__ for more information). 
}
\arguments{
  \item{habitatDF}{
data.frame of habitat units (cells) across the landscape, with X and Y coordinates, CellID, and habitat covariates.
}
  \item{CellID}{
Vector of CellIDs for animal relocations, in sequential order, with 
NA values when the animal was missing.
}
  \item{rsfCovars}{
Character string of covariate names for modeling habitat selection. Include 'distance' if step length is to modeled simultenously with the habitat selection coefficients (iSSF). Otherwise, the discrete choice, or 
conditional logistic regression model is fit, and probability of selection isn't influenced by distance 
required to travel. 
}
  \item{probDetCovars}{
Character string of covariate names to be used to model the probability of detection. 
Set to NULL for fitting a standard iSSF.
}
  \item{distColumns}{
Names of X and Y coordinates in habitat dataframe, or NULL if not including distance from previous location (aka step length).
}
  \item{maxLagArg}{
Recommended that you Do Not Disturb. Right now the cutoff is 4 missing observations in a row. If there are 5 or more, then the algorithm pretends that there are only 4. This is done to save computation time, and maybe improve convergence. However, it should be checked how often >4 consecutive missing fixes occur. Maybe play with this arguement a bit to see what could be gained. Likely very little if there are few really long strings of missing fixes in your data. 
}
}
\value{
Results of model estimation.
}
\seealso{
\code{\link{PDiSSF}}
\code{\link{habitat}}
\code{\link{locations}}
}
\examples{
  \dontrun{
  data(habitat)
  data(locations)
  
  # Fix success rate
  mean(!is.na(locations$unitID))
  
  # Computing time for larger data sets may vary WIDELY
  
  # Standard conditional logistic, if distColumn = NULL.
  pdissf(habitatDF = habitat, 
      CellID = locations$unitID, 
      iSSFCovars = c("prctSage", "elevation"), 
      probDetCovars = NULL)
  
  # integrated step selection function (iSSF) using step length,
  pdissf(habitatDF = habitat, 
      CellID = locations$unitID, 
      iSSFCovars = c("distance", "prctSage", "elevation"), 
      probDetCovars = NULL, 
      distColumns = c("utmX","utmY"))
  
  # PDRSF including step length
  pdissf(habitatDF = habitat, 
      CellID = locations$unitID, 
      iSSFCovars = c("distance", "prctSage", "elevation"), 
      probDetCovars = "prctSage", distColumns = c("utmX","utmY"))
  }
}
