est.utpn <- function(y, x=NULL, type=1, link="logit", q=0.5){
if(any(y<=0 | y>=1)) stop("y must be between 0 and 1")
if(!any(link == c("logit"))) 
            stop("link is not recognized")
if(type!=1 & type!=2) stop("type must be 1 or 2")
if(is.null(x)) {x <- matrix(1, ncol=1, nrow=length(y)); colnames(x)<-"Intercept"}
r=ncol(x)
llike.utpn <- function(theta, y, x, type=1, link="logit", q=0.5){
beta=theta[1:r]; lambda=theta[r+1]
dist <- switch(link, logit="logis")
g <- get(paste("d", dist, sep = ""), mode = "function")
G <- get(paste("p", dist, sep = ""), mode = "function")
rho=G(x %*% beta)
if(type==1)
{
	log.sigma <- log1p(-rho)-log(rho)-log(lambda+qnorm(1-q*pnorm(lambda)))
	ll<--log.sigma-2*log(y)-pnorm(lambda,log.p=TRUE)+dnorm((1-y)*exp(-log.sigma)/y-lambda,log=TRUE)
}
if(type==2)
{
	log.sigma <- log1p(rho)-log1p(-rho)-log(lambda+qnorm(1-(1-q)*pnorm(lambda)))
	ll <- -log.sigma-2*log1p(-y)-pnorm(lambda,log.p=TRUE)+dnorm(y*exp(-log.sigma)/(1-y)-lambda,log=TRUE)
}
-sum(ll)}
aux <- optim(rep(0,r+1), llike.utpn, y=y, x=x, type=type, link=link, method="BFGS",
control = list(maxit = 10000))
param = matrix(aux$par, ncol=1)
colnames(param) = c("estimate")
se = try(skewMLRM::solve2(pracma::hessian(llike.utpn, x0 = param, y = y, x=x, type=type, link=link)), silent = TRUE)
llike = -aux$value
if (!grepl("Error", se)[1]) {
        if (min(diag(se)) > 0) {
            param = cbind(param, sqrt(diag(se)))
            colnames(param) <- c("estimate", "s.e.")
        }
    }
rownames(param) <- c(paste("beta",colnames(x),sep=""), "lambda")
res = list(estimate = param, dist=paste("utpn",type,sep=""), conv = aux$conv, logLik = llike, 
AIC = -2 * llike + 2 * 2, BIC = -2 * llike + 2 * log(length(y)))
if (ncol(param) == 1) res$warnings = "Standard errors can't be estimated"
res
}