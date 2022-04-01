est.btpn <-
function (y) 
{
llike.BTPN<-function(theta, y, trans.param = FALSE){
sigma=theta[1]
lambda=theta[2]
eta=theta[3]
if (trans.param) sigma = exp(theta[1])
ll=-log(2)-log(sigma)-pnorm(lambda, log.p=TRUE)+dnorm(y*sqrt(1+eta^2)/(sigma*(sqrt(1+eta^2)-eta))-lambda,log=TRUE)*ifelse(y>=0,1,0)+dnorm(-y*sqrt(1+eta^2)/(sigma*(sqrt(1+eta^2)+eta))-lambda,log=TRUE)*ifelse(y<0,1,0)
-sum(ll)
}
aux = optim(c(0, 0, 0), llike.BTPN, y = y, trans.param = TRUE, 
        method = "Nelder-Mead", control = list(maxit = 10000))
    param = cbind(c(exp(aux$par[1]), aux$par[2:3]))
    colnames(param) = c("estimate")
    se = try(solve2(hessian(llike.BTPN, x0 = param, y = y)), silent = TRUE)
    llike = -aux$value
    if (!grepl("Error", se)[1]) {
        if (min(diag(se)) > 0) {
            param = cbind(param, sqrt(diag(se)))
            colnames(param) <- c("estimate", "s.e.")
        }
    }
    rownames(param) <- c("sigma", "lambda", "eta")
    res = list(estimate = param, conv = aux$conv, logLik = llike, 
        AIC = -2 * llike + 2 * 3, BIC = -2 * llike + 3 * log(length(y)))
    if (ncol(param) == 1) 
        res$warnings = "Standard errors can't be estimated"
    res

}
