matrixCalcFastTwo <- function (matrixOne, multTimes, matrixOfDiag, numberOfMultiplications) {
  .Call("_PDiSSF_matrixCalcFastTwo", PACKAGE = "PDiSSF", matrixOne, 
        multTimes, matrixOfDiag,numberOfMultiplications)
}
