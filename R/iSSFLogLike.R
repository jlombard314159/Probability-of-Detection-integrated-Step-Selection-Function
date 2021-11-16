#' @name iSSFLogLike
#'
#' @description The log-likelihood formulation for the iSSF
#' 
#' @param beta Coefficients to be estimated
#' 
#' @param X1 Data from the covariates for the iSSF component
#' 
#' @param locations data frame of GPS fixes
#' 
#' @param k1 Number of covariates for the iSSF
#' 


iSSFLogLike <- function (beta, X1, locations, k1) {
  
  iSSF.coefs <- beta[1:k1]
  ncells <- nrow(X1)
  fix.attempts <- length(locations)
  Lag <- rep(1, times = fix.attempts)
  for (i in 3:fix.attempts) {
    if (is.na(locations[i - 1])) 
      Lag[i] <- Lag[i - 1] + 1
  }
  detected <- !is.na(locations)
  EXP <- matrix(0, nrow = ncells, ncol = ncells)
  for (i in 1:k1) {
    EXP <- EXP + X1[, (1 + (i - 1) * ncells):(i * ncells)] * 
      iSSF.coefs[i]
  }
  EXP <- exp(EXP)
  sumEXP <- matrix(rowSums(EXP), nrow = ncells, ncol = ncells)
  D <- (EXP/sumEXP)
  prob <- rep(NA, times = fix.attempts)
  for (i in 2:fix.attempts) {
    if (detected[i] & Lag[i] == 1) {
      prob[i] <- D[locations[i - 1], locations[i]]
    }
  }
  cat(paste("Progress made at", Sys.time(), "\r"))
  -sum(log(prob), na.rm = TRUE)
}