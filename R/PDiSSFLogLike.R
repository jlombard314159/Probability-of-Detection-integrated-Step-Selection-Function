PDiSSFLogLike <- function (beta, X1, X2, locations, k1, k2, maxLagArg) {
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
        if (length(3:j) > maxLagArg) {
          numMults <- maxLagArg
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
  print(paste("meaningless", round(runif(1, 0, 100), 0)))
  -sum(log(prob), na.rm = TRUE)
}
