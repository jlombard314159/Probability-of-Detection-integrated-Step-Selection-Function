#' @name pdissf
#'
#' @description Helper function to set up the log-likelihood estimation
#' 
#' @param habitatDF Dataframe of habitat information
#' 
#' @param CellID Column that contains the unique ID for a cell location in the study
#' area
#' 
#' @param habitatCellID Column in the habitat DF that contains unique ID for the
#' habitat data
#' 
#' @param maxLagArg Number of matrix multiplications to toggle the matrix by
#' squaring algorithm. The default is 4 and it is not recommended to change this.
#' 
#' @param iSSFCovars covariates in the integrated step selection function model
#' 
#' @param probDetCovars covariates for the probability of detection
#' 
#' @param distColumns If distance is needed as a covariate, in either iSSF or the
#' probability of detection, then identify the columns that contain this information
#' 

pdissf <- function(habitatDF, CellID,  
                   maxLagArg, habitatCellID = 'unitID',
                   iSSFCovars = NULL, 
                 probDetCovars = NULL, distColumns = NULL) {
  
  distMatrix <- NULL
  
  nCells <- nrow(habitatDF)
  if (!is.data.frame(habitatDF)) {
    stop(paste("habitatDF should be a data frame."))
  }
  
  if(!(habitatCellID %in% colnames(habitatDF))){
    
    stop(paste('Your habitat DF does not have the column: ', habitatCellID))
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
                maxLagArg = maxLagArg, iSSFCovar = iSSFCovars, 
                LogCovar = c("Intercept", probDetCovars))
  return(PDiSSFFit)
}
