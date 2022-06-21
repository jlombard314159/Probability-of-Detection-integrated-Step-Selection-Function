#' pdissf
#' 
#' @name pdissf
#' 
#' @title Fits the probability of detection integrated step selection function (PDiSSF)
#'
#' @description This package fits the probability of detection integrated detection 
#' function presented in Vales et al. (2022) and presented as the "RSF for GPS fix success" 
#' in Nielson et al. (2009). As of now, this package relies on 
#' totally pre-processed data. See help(habitat) and help(locations) for more 
#' details and examples.
#' 
#' @param habitatDF Dataframe of habitat information
#' 
#' @param habitatCellID Column in the habitat DF that contains unique IDs for the 
#' habitat data. Some will have been used and thus also in the vector of CellID. 
#' Missing values are not allowed.
#' 
#' @param CellID Column name in the locations file that contains the ID for cells 
#' in the study area as they were used by the animal (chronological order is important).
#'  Should have IDs found in habitatCellID. Missing fixes should be 'NA'.
#' 
#' @param maximumGap Maximum allowable number of consecutive missing fixes. 
#' Default is 3 (3 consecutive missing locations). Large gaps will cause an 
#' exponential increase in memory requirements and computing time. This may 
#' result in convergence issues.
#' 
#' @param iSSFCovars covariates for step selection
#' 
#' @param probDetCovars covariates for the probability of detection
#' 
#' @param distColumns If distance from previous location (step length) is to be 
#' included in the model, identify the columns in the habitat DF that contain x-y 
#' coordinates in 2-D planar coordinates such as UTM or State Plane.
#' 
#' 
#' @usage pdissf(
#' habitatDF,
#' habitatCellID,
#' CellID,
#' maximumGap = 3,
#' iSSFCovars = NULL,
#' probDetCovars = NULL,
#' distColumns = NULL
#' )
#' 
#' @examples
#' library(PDiSSF)
#' data(habitat)
#' data(locations)
#' 
#' # Fix success rate
#' mean(!is.na(locations$unitID))
#' 
#' # Computing time for larger data sets may vary
#' 
#'# integrated step selection function (iSSF) using step length,
#'pdissf(habitatDF = habitat,
#'       CellID = locations$unitID, 
#'       habitatCellID = 'unitID',
#'       iSSFCovars = c('distance', 'prctSage', 'elevation'), 
#'       probDetCovars = NULL, 
#'       distColumns = c('utmX','utmY'))
#'
#'# PDiSSF including step length
#'pdissf(habitatDF = habitat, 
#'       CellID = locations$unitID, 
#'       habitatCellID = 'unitID',
#'       iSSFCovars = c('distance', 'prctSage', 'elevation'), 
#'       probDetCovars = 'prctSage', distColumns = c('utmX','utmY'),
#'       maximumGap = 3)
#' 
#' 
pdissf <- function(habitatDF, CellID,
                   habitatCellID,
                   maximumGap = 3, 
                   iSSFCovars = NULL, 
                   probDetCovars = NULL, 
                   distColumns = NULL) {
  
  distMatrix <- NULL
  
  nCells <- nrow(habitatDF)
  if (!is.data.frame(habitatDF)) {
    stop(paste("habitatDF should be a data frame."))
  }
  
  if(!(habitatCellID %in% colnames(habitatDF))){
    
    stop(paste('Your habitat DF does not have the column: ', habitatCellID))
  }
  
  if(maximumGap <1){
    stop(paste('Maximum Gap can not be less than 1. It must be a positive integer.'))
  }
  
  noNACellID <- unique(CellID[!is.na(CellID)])
  
  locationNotInHabitat <- noNACellID %in% habitatDF[,habitatCellID]
  
  idMismatch <- any(!locationNotInHabitat)
  
  if(idMismatch){
    errorMessage <- paste("at least one CellID does not match any of the habitatIDs. Check the following IDs: ",
                          noNACellID[!locationNotInHabitat])
    
    stop(errorMessage)
    
  }

  #Account for non-numeric cell ID
  habitatDF <- convertUnitIDToNumeric(habitatDF, unitCol = habitatCellID)
  
  updatedLocations <- updateLocationUnitID(habitatData=habitatDF,
                                            locationData = CellID,
                                            habitatUnitCol = habitatCellID)
  
  CellID <- updatedLocations$numericUnitID
  
  habitatDF[,habitatCellID] <- habitatDF$numericUnitID

  factorHabitat <- sapply(habitatDF, is.factor)
  factorHabitat <- factorHabitat[!(factorHabitat == FALSE)]
  
  toConvertHabitat <- names(factorHabitat)
  factorDF <- as.data.frame(habitatDF[, colnames(habitatDF) %in% 
                                        toConvertHabitat])

  habitatDF <- habitatDF[,sort(names(habitatDF))] # just added
  iSSFCovars <- sort(iSSFCovars) # just added
  if(sum(iSSFCovars == "distance") == 1){ # just added to stick distance at the end
    iSSFCovars <- c(iSSFCovars[iSSFCovars != "distance"], "distance")
  }

  if (ncol(factorDF) > 0) {
    factorDFList <- list()
    for (i in 1:ncol(factorDF)) {
      mm <- model.matrix(~factorDF[, i] - 1, model.frame(~factorDF[, 
                                                                   i] - 1), contrasts = FALSE)
      factorDFList[[i]] <- cbind(mm[, 2:ncol(mm)])
      colnames(factorDFList[[i]]) <- gsub(".*]", "", colnames(factorDFList[[i]]))
    }
    newDF <- do.call(cbind.data.frame, factorDFList)
    if (toConvertHabitat %in% iSSFCovars) {
      iSSFCovars <- iSSFCovars[!iSSFCovars == toConvertHabitat]
      iSSFCovars <- c(iSSFCovars, colnames(newDF))
    }
    if (toConvertHabitat %in% probDetCovars) {
      probDetCovars <- probDetCovars[!probDetCovars == 
                                       toConvertHabitat]
      probDetCovars <- c(probDetCovars, colnames(newDF))
    }
    habitatDF <- habitatDF[, !(colnames(habitatDF) %in% toConvertHabitat)]
    habitatDF <- cbind(habitatDF, newDF)
  }
  
  if(!is.null(distColumns)){
  distMatrix <- as.matrix(x = stats::dist(cbind(habitatDF[,colnames(habitatDF) %in% distColumns[1]],
                                                habitatDF[,colnames(habitatDF) %in% distColumns[2]])))
  
  }
  
  if (is.null(iSSFCovars)) {
    stop("A covariate must be specific for the iSSF")
  } else {
    habitatDFSub <- as.data.frame(habitatDF[, colnames(habitatDF) %in% 
                                              iSSFCovars])
    matrixList <- list()
    for (i in 1:ncol(habitatDFSub)) {
      matrixList[[i]] <- matrix(habitatDFSub[, i], nrow = nCells, 
                                ncol = nCells, byrow = T)
    }

    if (!is.null(distMatrix) & ("distance" %in% iSSFCovars)) {
      matrixList[[length(matrixList) + 1]] <- distMatrix
    }
    vectorToConvert <- c()
    for (i in 1:length(matrixList)) {
      vectorToConvert[i] <- paste("matrixList[[", i, "]]", 
                                  sep = "")
    }
    formulaRest <- paste("~", paste(vectorToConvert, collapse = "+"))
    selectionFormula <- as.formula(formulaRest)
  }
  if (is.null(probDetCovars)) {
    probDetFormula <- NULL
  } else {
    logisticSub <- as.data.frame(habitatDF[, colnames(habitatDF) %in% 
                                             probDetCovars])
    vectorToConvert <- c()
    for (i in 1:ncol(logisticSub)) {
      vectorToConvert[i] <- paste("logisticSub[,", i, "]", 
                                  sep = "")
    }
    finalFormula <- paste("~ 1 +", paste(vectorToConvert, 
                                         collapse = "+"))
    probDetFormula <- as.formula(finalFormula)
  }
  PDiSSFFit <- PDiSSFModel(selection = selectionFormula, p = probDetFormula, 
                locations = CellID, ncells = nCells, 
                maximumGap = maximumGap, iSSFCovar = iSSFCovars, 
                LogCovar = c("Intercept", probDetCovars))
  return(PDiSSFFit)
}
