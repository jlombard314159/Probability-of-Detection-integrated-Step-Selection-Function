convertUnitIDToNumeric <- function(habitatData, unitCol='unitID'){
  
  sortedVector <- sort(habitatData[,unitCol])
  
  indexVector <- seq(1,length(sortedVector),1)
  
  mergeDF <- data.frame(unitCol = sortedVector, 
                        'numericUnitID'=indexVector)
  
  habitatData$numericUnitID <- mergeDF$numericUnitID[match(paste(habitatData$unitID),
                                                           paste(mergeDF$unitCol))]
  
  
  return(habitatData) 
  
}

updateLocationUnitID <- function(habitatData, locationData, habitatUnitCol = 'unitID'){
  
  locationDataFrame <- data.frame('unitID' = locationData)
  
  locationDataFrame[,'numericUnitID'] <- habitatData$numericUnitID[match(paste(locationDataFrame$unitID),
                                                                         paste(habitatData[,habitatUnitCol]))]
  
  
  return(locationDataFrame)
  
}