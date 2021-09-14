PDiSSFModel <- function (selection, p, locations, ncells, maxLagArg, iSSFCovar, 
    LogCovar,numberOfMatrixMult) {
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
                (k.iSSF), bic = 2 * out$objective + (k.iSSF) * 
                log(sum(!is.na(locations)) - 1))
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
            k2 = k.p, maxLagArg = maxLagArg,
            numberOfMatrixMult = numberOfMatrixMult)
        
        iSSF.coefs <- out$par[1:k.iSSF]
        p.coefs <- out$par[(k.iSSF + 1):(k.iSSF + k.p)]
        
        hessian <- F.2nd.deriv(out$par, PDiSSFLogLike, X1 = X.iSSF, 
                    X2 = X.p, locations = locations, k1 = k.iSSF, 
                    k2 = k.p, maxLagArg = maxLagArg,
                    numberOfMatrixMult = numberOfMatrixMult)
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
