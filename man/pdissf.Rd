\name{pdissf}
\alias{pdissf}
\title{
Fits the probability of detection integrated step selection function (PDiSSF)
}

\description{

This package fits the probability of detection integrated detection function presented in Vales et al. (2022) and presented as the "RSF for GPS fix success" in Nielson et al. (2009). As of now, this package relies on 
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
data.frame of habitat units (cells) across the analysis area (e.g., home range), with X and Y coordinates, CellID, and habitat covariates. The X and Y coordinates need to be in UTM. 
}
  \item{CellID}{
Vector of CellIDs for animal relocations, in sequential order of use, with 
NA values when the GPS unit didn't record a scheduled fix. The sequential location data must start and end with successful fixes. To emphasize, order your location data by date and time first. 

CellID values must match those in the habitatDF.
}
  \item{issfCovars}{
Character string of covariate names for modeling habitat selection. Include 'distance' in the list of habitat selection coefficients if you want to include step length in the model. If 'distance' is not included, the discrete choice, or conditional logistic regression, model is fit considering the entire study area to be equally available, and the probability of selection is not influenced by the distance required to travel from location A to B. 
}
  \item{probDetCovars}{
Character string of covariate names for estimating the probability of detection. 
Set to NULL for fitting a standard iSSF.
}
  \item{distColumns}{
Names of X and Y coordinates in habitat data, or left as NULL if not including distance from the previous location (aka step length).
}
  \item{maxLagArg}{
When there is a large gap of consecutive missing locations (Lag) in the data, model fitting can be slow, and you may experience convergence issues. By default, the algorithm is trying to find where the animal could have gone between GPS fixes i and i + (Lag + 1), considering all possible paths between the two successful locations bounding the gap in the data. The maxLagArg can be used to limit the number of possible locations between locations i and i + (Lag + 1). For example, suppose there are 7 missing locations in a row and maxLagArg = 5. In that case, the algorithm finds the most likely path between the two successful locations pretending there are only 5 missing locations instead of 7. Use this option with caution, as setting a maxLagArg < the maximum lag in the location data will result in positively biased estimates of the average step length (distance traveled). The default is the maximum lag in the data, but maxLagArg can be modified for testing purposes. 
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
As-is, the analysis ignores any changes in fix schedule. 

The computing demand can be high for larger data sets, depending on fix success, the number of available cells, the number of fix attempts, and the number of covariates in the model. Fitting the PDiSSF could take seconds or days. The size of the study area (number of cells) has the most influence on memory requirements, processing time, and convergence. Parallel processing (try library parallel) can help estimate models for multiple animals; multithreading with Microsoft R Open can speed up processing. However, parallel processing will not speed up model estimation but merely allow the user to fit more than one model at a time. See the Supporting Information in Vales et al. (2022; __github__ for more information). 

The system time is printed to the console to indicate that model estimation is ongoing and can be used to identify when the process has stalled. The system time is not printed when running in parallel.
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
  # calculate maximum number consecutive missing locations 
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
