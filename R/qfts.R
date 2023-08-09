qfts <-
function(p, sigma=1, lambda=1, dist="norm")
{
if(!any(dist == c("norm", "logis", "cauchy", "laplace"))) 
            stop("distribution is not recognized")
if(any(p<0 | p>1)) stop("p must be between zero and one")
if(sigma<=0) stop("sigma must be positive")
G <- get(paste("p", dist, sep = ""), mode = "function")
Q <- get(paste("q", dist, sep = ""), mode = "function")
sigma*(lambda+Q(1-(1-p)*G(lambda)))
}
