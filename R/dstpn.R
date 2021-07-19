dstpn <-
function(x, sigma, lambda, q, log=FALSE)
{
if(is.null(x))
    stop("x must be specified")
  if(is.null(sigma))
    stop("sigma must be specified")
  if(is.null(lambda))
    stop("lambda must be specified")
  if(is.null(q))
    stop("q must be specified")
  if(sigma<=0)
    stop("sigma must be positive")
  if(q<=0)
    stop("q must be positive")
  if(any(x<=0))
    stop("x's must be positive")
aux.int=function(x, y, lambda, sigma, q) exp(q*log(x)+dnorm(x*y/sigma-lambda,log=TRUE))
int=c(); for(i in 1:length(x)){
int[i]=integrate(aux.int, lower=0, upper=1, y=x[i], lambda=lambda, sigma=sigma, q=q)$value}
lf<-log(q)-log(sigma)-pnorm(lambda, log.p=TRUE)+log(int)
  if(!log) lf=exp(lf)
  lf
}
