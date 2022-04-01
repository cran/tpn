dftpn <-
function(x, sigma=1, lambda=1, dist="norm", log=FALSE)
{
if(!any(dist == c("norm", "logis", "cauchy", "laplace"))) 
            stop("distribution is not recognized")
#if(any(x<=0)) stop("x must be positive")
if(sigma<=0) stop("sigma must be positive")
g <- get(paste("d", dist, sep = ""), mode = "function")
G <- get(paste("p", dist, sep = ""), mode = "function")
lf <- -log(sigma)-G(lambda, log.p=TRUE)+g(x/sigma-lambda,log=TRUE)
if(!log) lf <- exp(lf)
lf
}
