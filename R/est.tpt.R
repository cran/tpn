est.tpt=function (y, x = NULL, q = 0.5) 
{
    if (any(y <= 0) ) 
        stop("y must be positive")
    if (is.null(x)) {
        x <- matrix(1, ncol = 1, nrow = length(y))
        colnames(x) <- "Intercept"
    }
    if (!is.matrix(x)) 
        stop("x must be a matrix")
    if (nrow(x) != length(y)) 
        stop("rows of x don't agree with the length of y")
    r = ncol(x)
    if (is.null(colnames(x))) 
        colnames(x) <- paste("var", 1:r, sep = "")
    if (r > 1) {
        llike.tpt <- function(theta, y, x, q = 0.5, t.param = FALSE) {
            beta = theta[1:r]
            lambda = theta[r + 1]
		nu=theta[r+2]
		if(t.param) nu=exp(theta[r+2])
            rho = exp(x %*% beta)
		sigma = rho/(lambda+qt(1+pt(lambda, df=nu)*(q-1),df=nu))
		ll = -log(sigma) - pt(lambda, df = nu, log.p = TRUE) + 
            dt(y/sigma - lambda, df = nu, log = TRUE)
            -sum(ll)
        }
        aux <- optim(rep(0, r + 2), llike.tpt, y = y, x = x, 
            q = q, method = "BFGS", control = list(maxit = 10000), t.param=TRUE)
        param = matrix(c(aux$par[-(r+2)],exp(aux$par[r+2])), ncol = 1)
        colnames(param) = c("estimate")
        llike = -aux$value
	  se = try(solve2(hessian(llike.tpt, x0 = param, y = y, x=x)), silent = TRUE)
        if (!grepl("Error", se)[1]) {
            if (min(diag(se)) > 0) {
                param = cbind(param, sqrt(diag(se)))
                colnames(param) <- c("estimate", "s.e.")
            }
        }
        rownames(param) <- c(paste("beta", colnames(x), sep = ""), 
            "lambda","nu")
        res = list(estimate = param, dist = "tpt", quantile = q, 
		conv = aux$conv, logLik = llike, 
            AIC = -2 * llike + 2 * nrow(param), BIC = -2 * llike + 
                nrow(param) * log(length(y)))
        if (ncol(param) == 1) 
            res$warnings = "Standard errors can't be estimated"
    }
    if (r == 1) {
        llike.tpt0 <- function(theta, y, t.param = FALSE) {
        sigma = theta[1]
        lambda = theta[2]
        nu = theta[3]
        if (t.param) {
            sigma = exp(theta[1])
            nu = exp(theta[3])
        }
        ll = -log(sigma) - pt(lambda, df = nu, log.p = TRUE) + 
            dt(y/sigma - lambda, df = nu, log = TRUE)
        -sum(ll)
    }
    aux = optim(c(0, 0, 0), llike.tpt0, y = y, t.param = TRUE, 
        method = "BFGS", control = list(maxit = 10000))
    param = cbind(c(exp(aux$par[1]), aux$par[2], exp(aux$par[3])))
    colnames(param) = c("estimate")
    se = try(solve2(hessian(llike.tpt0, x0 = param, y = y)), silent = TRUE)
    llike = -aux$value
    if (!grepl("Error", se)[1]) {
        if (min(diag(se)) > 0) {
            param = cbind(param, sqrt(diag(se)))
            colnames(param) <- c("estimate", "s.e.")
        }
    }
    rownames(param) <- c("sigma", "lambda", "nu")
    res = list(estimate = param, conv = aux$conv, logLik = llike, 
        AIC = -2 * llike + 2 * 3, BIC = -2 * llike + 3 * log(length(y)))
    if (ncol(param) == 1) 
        res$warnings = "Standard errors can't be estimated"
    }
    res
}