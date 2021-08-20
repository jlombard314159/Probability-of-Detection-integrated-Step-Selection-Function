matrixCalcFastTwo <- function (matrixOne, multTimes, matrixOfDiag) {
  .Call("_PDiSSF_matrixCalcFastTwo", PACKAGE = "PDiSSF", matrixOne, 
        multTimes, matrixOfDiag)
}
