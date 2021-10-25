\name{pdissf}
\alias{pdissf}
\title{
Fits the probability of detection integrated step selection function (PDiSSF)
}

\description{

This package fits the probability of detection integrated detection function presented in Vales et al. (2022), and presented as the "RSF for GPS fix success" in Nielson et al. (2009). As of now, this package relies on 
totally pre-processed data. See help(habitat) and help(locations) for more details and examples.
}

\usage{
pdissf(habitatDF = habitat, 
      CellID = locations$unitID, 
      iSSFCovars = c("distance", "prctSage", "elevation"), 
      probDetCovars = "prctSage", 
      distColumns = c("utmX","utmY"), 
      maxLagArg = 10) # needs to be updated
}

\arguments{
  \item{habitatDF}{
data.frame of habitat units (cells) across the analysis area (e.g., home range), with X and Y coordinates, CellID, and habitat covariates.
}
  \item{CellID}{
Vector of CellIDs for animal relocations, in sequential order, with 
NA values when the GPS unit was missing a scheduled fix. The sequential location data must start and end with successful fixes. 
}
  \item{issfCovars}{
Character string of covariate names for modeling habitat selection. Include 'distance' if step length is to be modeled simultaneously with the habitat selection coefficients (iSSF). If 'distance' is not included, the discrete choice, or 
conditional logistic regression model is fit, and the modeled probability of selection is not influenced by distance 
required to travel from location A to B. 
}
  \item{probDetCovars}{
Character string of covariate names to be used to model the probability of detection. 
Set to NULL for fitting a standard iSSF.
}
  \item{distColumns}{
Names of X and Y coordinates in habitat dataframe, or NULL if not including distance from previous location (aka step length).
}
  \item{maxLagArg}{
WE ARE REWORKING THIS
Right now the cutoff is 4 consecutive missing observations. If there are 5 or more, then the algorithm pretends that there are only 4. This is done to save computation time, and maybe improve convergence. However, it should be checked how often >4 consecutive missing fixes occur. Maybe play with this arguement a bit to see what could be gained. Likely very little if there are few really long strings of missing fixes in your data. 
}
}
\value{
Results of model estimation. 
Show results here?
convergence: 
  0 = converged
  1 = no convergence
  
}
\details{
The computing demand can be high for larger data sets, depending fix success, the number of available cells, the number of fix attempts, and the number of covariates in the model. Fitting the PDiSSF could take seconds or days. The size of the study area (number of cells) has the most influence on memory requirements, processing time, and convergence. Parallel processing (try library parallel) can help when analyzing multiple animals, multitheading with Mircrosf R Open can speed up procesing. Parallel processing will not speed up model estimation, but merely allow the user to fit more than one model at a time. See the Supporting Information in Vales et al. (2022; __github__ for more information). 
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
  
  # integrated step selection function (iSSF) using step length
  pdissf(habitatDF = habitat, 
      CellID = locations$unitID, 
      iSSFCovars = c("distance", "prctSage", "elevation"), 
      probDetCovars = NULL, 
      distColumns = c("utmX","utmY"))
  
  # problem with hessian: center elevation to fix problem
  habitat$elevationCntrd <- habitat$elevation - mean(habitat$elevation)
  pdissf(habitatDF = habitat, 
      CellID = locations$unitID, 
      iSSFCovars = c("distance", "prctSage", "elevationCntrd"), 
      probDetCovars = NULL, 
      distColumns = c("utmX","utmY"))

  # PDiSSF with prctSage as covariate for probability of detection
  # calculate maximum number consecutive missing locations MODIFYING THIS
  missing <- ifelse(is.na(locations$unitID), 1, 0)
  nAttempts <- length(locations$unitID)
  sinceFix <- rep(0, nAttempts)
  for(i in 2:nAttempts){
    if(missing[i] == 1){
      sinceFix[i] <- missing[i] + sinceFix[i-1] 
    } else {
      sinceFix[i] <- 0
    }
  }
  max(sinceFix) # set maxLagArg to 3 or lower
  
  pdissf(habitatDF = habitat, 
      CellID = locations$unitID, 
      iSSFCovars = c("distance", "prctSage", "elevation"), 
      probDetCovars = "prctSage", distColumns = c("utmX","utmY"), maxLagArg = 3)

  
  }
}