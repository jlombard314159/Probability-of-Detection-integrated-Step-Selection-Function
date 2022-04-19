#' PDiSSFLogLike
#' 
#' @name PDiSSFLogLike
#'
#' @description Calculates the loglikelihood value for a set of parameters
#' 
#' @param beta Coefficients to be estimated
#' 
#' @param X1 Data from the covariates for the iSSF component
#' 
#' @param X2 Data from the covariates for the probability of detection component
#' 
#' @param locations data frame of GPS fixes
#' 
#' @param k1 Number of covariates for the iSSF
#' 
#' @param k2 Number of covariates for probability of detection. Equal to 0
#' if no probability of detection is being estimated.
#' 
#' @param maximumGap Maximum allowable number of consecutive missing fixes. Default is 4
#' (3 consecutive missing locations). Larger max lags will increase computing time
#' and may result in non convergence.
#' 

PDiSSFLogLike <- function (beta, X1, X2, locations, k1, k2, maximumGap) {
  link <- function(x) {
    exp(x)/(1 + exp(x))
  }
  iSSF.coefs <- beta[1:k1]
  p.coefs <- t(t(beta[(k1 + 1):(k1 + k2)]))
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
  P <- link(X2 %*% p.coefs)
  D <- (EXP/sumEXP) * matrix(P, nrow = ncells, ncol = ncells, 
                            byrow = TRUE)
  B <- (EXP/sumEXP) * matrix(1 - P, nrow = ncells, ncol = ncells, 
                            byrow = TRUE)
  if (max(Lag) > 1) {
    B.s <- vector("list", max(Lag - 1))
    B.s[[1]] <- diag(ncells)
    B.s[[2]] <- B
    if (max(Lag) > 2) {
      for (j in 3:max(Lag - 1)) {
        output <- B
        if (length(3:j) > maximumGap) {
          numMults <- maximumGap
        }
        else {
          numMults <- length(3:j)
        }
        
        B.s[[j]] <- matrixCalcFastTwo(matrixOne = output,
                                      multTimes = numMults, 
                                      matrixOfDiag = diag(1, ncells))
      }
    }
  }
  prob <- rep(NA, times = fix.attempts)
  for (i in 2:fix.attempts) {
    if (detected[i] & Lag[i] == 1) {
      prob[i] <- D[locations[i - 1], locations[i]]
    }
    if (detected[i] & Lag[i] > 1) {
      prob[i] <- B[locations[i - Lag[i]], ] %*% B.s[[Lag[i] - 
                                                      1]] %*% D[, locations[i]]
    }
  }
  cat(paste("Progress made at", Sys.time(), "\r"))
  -sum(log(prob), na.rm = TRUE)
}
