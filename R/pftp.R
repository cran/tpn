pftp <-
function(x, sigma=1, lambda=1, dist="norm", lower.tail=TRUE, log.p=FALSE)
{
if(!any(dist == c("norm", "logis", "cauchy", "laplace"))) 
            stop("distribution is not recognized")
if(any(x<=0)) stop("x must be positive")
if(sigma<=0) stop("sigma must be positive")
g <- get(paste("d", dist, sep = ""), mode = "function")
G <- get(paste("p", dist, sep = ""), mode = "function")
lf <- G(x/sigma-lambda, log.p=TRUE)-G(lambda, log.p=TRUE, lower.tail=FALSE)-G(lambda, log.p=TRUE)
if(log.p) lf <- exp(lf)
if(!lower.tail) lf <- 1-lf
lf
}
