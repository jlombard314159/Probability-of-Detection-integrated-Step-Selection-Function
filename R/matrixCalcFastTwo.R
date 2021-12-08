#' matrixCalcFastTwo
#' 
#' @name matrixCalcFastTwo
#'
#' @description Calls the rcpp file
#' 
#' @param matrixOne Output to be multiplied multTimes times
#' 
#' @param multTimes number of matrix multiplications to perform
#' 
#' @param matrixOfDiag Empty matrix of the form diag(1, ncells)
#' 
#' 

matrixCalcFastTwo <- function (matrixOne, multTimes, matrixOfDiag) {
  .Call("_PDiSSF_matrixCalcFastTwo", PACKAGE = "PDiSSF", matrixOne, 
        multTimes, matrixOfDiag)
}
