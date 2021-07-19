rtpn <-
function(n, sigma, lambda)
{
if(is.null(n))
    stop("sample size must be specified")
  if(is.null(sigma))
    stop("sigma must be specified")
  if(is.null(lambda))
    stop("lambda must be specified")
  if(n<=0 | round(n)!=n)
    stop("sample size must be a positive integer")      
v=runif(n)
u=runif(n)
sigma*(qnorm(1+pnorm(lambda)*(v-1))+lambda)
}
