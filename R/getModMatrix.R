#' getModMatrix
#' 
#' @name getModMatrix
#'
#' @description Mainly a helper function, this creates the model matrix for
#' the optimization process. Not recommended for isolated use.
#' 
#' @param f model formula
#' 
#' @param ncells number of cells in the study area


getModMatrix <- function (f, ncells){
  call <- match.call()
  contrasts <- NULL
  mf <- match.call(expand.dots = FALSE)
  mf$family <- mf$start <- mf$control <- mf$maxit <- NULL
  mf$model <- mf$method <- mf$x <- mf$y <- mf$contrasts <- NULL
  mf$... <- NULL
  mf$f <- mf$ncells <- NULL
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf$formula <- f
  form <- as.character(formula(mf))[2]
  intercept <- NA
  if (nchar(form) == 1 & form == "1") {
    X <- matrix(1, ncells, 1)
    intercept <- TRUE
    nx <- 1
    x.names <- character(0)
  }
  else {
    if (substr(form, start = 1, stop = 1) == "1") {
      int <- 1
      intercept <- TRUE
    }
    else {
      int <- 0
    }
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    if (int == 0) 
      attr(mt, "intercept") <- 0
    xvars <- as.character(attr(mt, "variables"))[-1]
    if ((yvar <- attr(mt, "response")) > 0) 
      xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
      xlev <- lapply(mf[xvars], levels)
      xlev[!sapply(xlev, is.null)]
    }
    X <- NA
    X <- if (length(attr(mt, "order")) != 0) {
      model.matrix(mt, mf, contrasts)
    }
    else {
      stop("Empty models not allowed")
    }
    assign.col <- attr(X, "assign")
    nx <- length(unique(assign.col))
    x.names <- xvars
    dimnames(X) <- NULL
  }
  ans <- list(X = X, n.covars = nx, intercept = intercept, 
              vars = x.names)
  return(ans)
}
