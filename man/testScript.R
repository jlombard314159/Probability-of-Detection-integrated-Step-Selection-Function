## Not run: 
library(PDiSSF)
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

## End(Not run)
