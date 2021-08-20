F.2nd.deriv <- 
  function (x, FUN, ...) 
  {
    FUN <- match.fun(FUN)
    d <- length(x)
    hess <- matrix(0, nrow = d, ncol = d)
    eps <- 1e-06
    h <- ifelse(x == 0, eps^0.25, (eps^(0.25)) * x)
    for (i in 1:d) {
      ei <- rep(0, d)
      ei[i] <- 1
      hess[i, i] <- (-FUN(x + 2 * h * ei, ...) + 16 * FUN(x + 
                                                            h * ei, ...) - 30 * FUN(x, ...) + 16 * FUN(x - h * 
                                                                                                         ei, ...) - FUN(x - 2 * h * ei, ...))/(12 * h[i] * 
                                                                                                                                                 h[i])
      if ((i + 1) <= d) {
        for (j in (i + 1):d) {
          ej <- rep(0, d)
          ej[j] <- 1
          hess[i, j] <- (FUN(x + h * ei + h * ej, ...) - 
                           FUN(x + h * ei - h * ej, ...) - FUN(x - h * 
                                                                 ei + h * ej, ...) + FUN(x - h * ei - h * ej, 
                                                                                         ...))/(4 * h[i] * h[j])
          hess[j, i] <- hess[i, j]
        }
      }
    }
    hess
  }
