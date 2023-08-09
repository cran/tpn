est.fts <-
function(y, dist="norm"){
if(any(y<=0)) stop("y must be positive")
if(!any(dist == c("norm", "logis", "cauchy", "laplace"))) 
            stop("distribution is not recognized")
lambda.ini=ifelse(moments::skewness(y)<0,-1,1)
sigma.ini=sqrt(mean(y^2))
llike.ftpn <- function(theta, x, dist="norm", t.param=FALSE){
sigma=theta[1]; lambda=theta[2]
if(t.param)sigma=exp(theta[1])
g <- get(paste("d", dist, sep = ""), mode = "function")
G <- get(paste("p", dist, sep = ""), mode = "function")
ll <- -log(sigma)-G(lambda, log.p=TRUE)+g(x/sigma-lambda,log=TRUE)
-sum(ll)}
aux <- optim(c(log(sigma.ini),lambda.ini), llike.ftpn, x=y, dist=dist, t.param=TRUE, method="Nelder-Mead",
control = list(maxit = 10000))
aux <- optim(aux$par, llike.ftpn, x=y, dist=dist, t.param=TRUE, method="BFGS",
control = list(maxit = 10000))
param = cbind(c(exp(aux$par[1]), aux$par[2]))
colnames(param) = c("estimate")
se = try(skewMLRM::solve2(pracma::hessian(llike.ftpn, x0 = param, x = y, dist=dist)), silent = TRUE)
llike = -aux$value
if (!grepl("Error", se)[1]) {
        if (min(diag(se)) > 0) {
            param = cbind(param, sqrt(diag(se)))
            colnames(param) <- c("estimate", "s.e.")
        }
    }
rownames(param) <- c("sigma", "lambda")
res = list(estimate = param, dist=dist, conv = aux$conv, logLik = llike, 
AIC = -2 * llike + 2 * 2, BIC = -2 * llike + 2 * log(length(y)))
if (ncol(param) == 1) res$warnings = "Standard errors can't be estimated"
res
}
