#' @docType data
#' @keywords exampleHabitat
#' @name exampleHabitat
#' @usage data(exampleHabitat)
#' @description A data.frame named habitat containing covariate values for all 
#' 30 m cells within the deer’s home range. Data from the analysis by Vales et al. (2022). 
#' See the manuscript for details. 
#' This data.frame contains the following information:
#'   \itemize{
#'    \item pointID: unique ID for each 30 m cell in the deer’s home range
#'    \item cc: Percent canopy cover for a cell
#'    \item cc75: cancopy cover shifted (centered) by subtracting 75 meters. 
#'    \item ab: accepted biomass (amount of neutral and selected forage)
#'    \item dcfe: distance (m) to cover-forage edge
#'    \item drds: distance (m) to nearest road
#'    \item ealbers: easting coordinate for the cell in UTM
#'    \item nalbers: northing coordinate for the cell in UTM
#'}

NULL

#' @docType data
#' @keywords locationHabitat
#' @name locationHabitat
#' @usage data(locationHabitat)
#' @description A data.frame named locations of GPS fixes from one deer used in 
#' the analysis by Vales et al. (2022). See the manuscript for details. 
#' This data.frame contains the following information:
#'    \itemize{
#'     \item lineno: Ignored. This is simply to recognize that all fix attempts 
#'        (missing and successful) need to be represented in the data. 
#'     \item The numbers should be consecutive and have no gaps.
#'     \item pointID: cell ID for 30 m cells where the deer was located at each fix attempt.
#'     \item Missing locations should be in the data and represented by NA. All pointID 
#'     \item values should be represented in the habitat data.frame. 
#'    }

NULL