rftp <-
function(n, sigma=1, lambda=1, dist="norm")
{
n<-round(n)
if(n<=0) stop("n must be positive")
if(!any(dist == c("norm", "logis", "cauchy", "laplace"))) 
            stop("distribution is not recognized")
if(sigma<=0) stop("sigma must be positive")
qftp(runif(n), sigma=sigma, lambda=lambda, dist=dist)
}
