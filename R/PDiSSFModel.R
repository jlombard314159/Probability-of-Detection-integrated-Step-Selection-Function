#' PDiSSFModel
#' 
#' @name PDiSSFModel
#'
#' @description Sets up the parameter estimation for the log-likelihood function.
#' Not recommended for isolated use.
#' 
#' @param selection Model for habitat selection
#' 
#' @param p Model for encounter probabilities
#' 
#' @param locations Dataframe of animal locations
#' 
#' @param ncells Number of cells in the study area
#' 
#' @param maximumGap Maximum allowable number of consecutive missing fixes. Default is 3
#' (3 consecutive missing locations). Larger max lags will increase computing time
#' and may result in non convergence.
#' 
#' @param iSSFCovar Covariates for the iSSF. Not explicitly specified by user, 
#' comes directly from pdissf.R.
#' 
#' @param LogCovar Covariates for the probability of detection, plus an intercept
#' term. Not explicitly specified by user, comes directly from pdissf.R
#' 
#' 

PDiSSFModel <- function (selection, p, locations, ncells, maximumGap, iSSFCovar, 
    LogCovar) {
    if (missing(selection)) 
        stop("Model for habitat selection must be specified")
    if (missing(p)) 
        stop("Model for p (encounter probabilities) must be specified")
    if (missing(locations)) 
        stop("Animal locations must be specified")
    if (max(locations, na.rm = T) > ncells) 
        stop("Cell ID #s are not correct")
    if (!is.numeric(locations)) # Consider moving one above
        stop("Cell ID #s are not numeric")
  
    orig.call <- match.call()
    iSSF.mod.mat <- getModMatrix(selection, ncells)
    X.iSSF <- iSSF.mod.mat$X
    k.iSSF <- iSSF.mod.mat$n.covars

    if (is.null(p)) {
        strt.vals <- rep(0, k.iSSF)
        out <- nlminb(start = strt.vals, objective = iSSFLogLike, 
            X1 = X.iSSF, locations = locations, k1 = k.iSSF)
        
        hessian = F.2nd.deriv(out$par, iSSFLogLike, X1 = X.iSSF, 
                         locations = locations, k1 = k.iSSF)
        
        SEs <- sqrt(diag(solve(hessian)))
        iSSF.coefs <- out$par[1:k.iSSF]
        iSSFDataframe <- data.frame(Covar = iSSFCovar, Coef = iSSF.coefs, 
            SE = SEs[1:k.iSSF])
        
        ans <- list(loglik = -out$objective, convergence = out$convergence, 
            call = orig.call, ncells = ncells, n.fix.attempts = length(locations), 
            iSSF = iSSFDataframe, aic = 2 * out$objective + 2 * 
                (k.iSSF))
        class(ans) <- "PDiSSF"
        return(ans)
    }
    if (!is.null(p)) {
        p.mod.mat <- getModMatrix(p, ncells)
        X.p <- p.mod.mat$X
        k.p <- p.mod.mat$n.covars
        strt.vals <- rep(0, k.iSSF + k.p)
        out <- nlminb(start = strt.vals, objective = PDiSSFLogLike, 
            X1 = X.iSSF, X2 = X.p, locations = locations, k1 = k.iSSF, 
            k2 = k.p, maximumGap = maximumGap)
        
        iSSF.coefs <- out$par[1:k.iSSF]
        p.coefs <- out$par[(k.iSSF + 1):(k.iSSF + k.p)]
        
        hessian <- F.2nd.deriv(out$par, PDiSSFLogLike, X1 = X.iSSF, 
                    X2 = X.p, locations = locations, k1 = k.iSSF, 
                    k2 = k.p, maximumGap = maximumGap)

        SEs <- sqrt(diag(solve(hessian)))

        # iSSFDataframe <- data.frame(Covar = iSSFCovar, Coef = iSSF.coefs[order(iSSFCovar)], ### original
        #     SE = SEs[1:k.iSSF])
        iSSFDataframe <- data.frame(Covar = iSSFCovar, Coef = iSSF.coefs, 
                                    SE = SEs[1:k.iSSF])
        LogDataframe <- data.frame(Covar = LogCovar, Coef = p.coefs, 
            SE = SEs[(k.iSSF + 1):(k.iSSF + k.p)])
        ans <- list(loglik = -out$objective, convergence = out$convergence, 
            call = orig.call, n.ells = ncells, n.fix.attempts = length(locations), 
            iSSF = iSSFDataframe, Detection = LogDataframe, aic = 2 * 
                out$objective + 2 * (k.iSSF + k.p))
        class(ans) <- "PDiSSF"
        return(ans)
    }
}
