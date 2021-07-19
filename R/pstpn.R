pstpn <-
function(x, sigma, lambda, q, lower.tail=TRUE, log=FALSE){
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
lf = function(x, sigma, lambda, q, xx) { exp(log(q)-log(sigma)-pnorm(lambda, log.p=TRUE)+q*log(x)+dnorm(x*xx/sigma-lambda,log=TRUE)) }
pf = function(y, sigma, lambda, q) { sapply(y,
    function(z) { integrate(lf, 0, 1, sigma=sigma, lambda=lambda, q=q, xx=z)$value }) }
pp=c();for(i in 1:length(x)){
pp[i]=integrate(pf, 0, x[i], sigma=sigma, lambda=lambda, q=q)$value}
if(!lower.tail) pp=1-pp
  if(log) pp=log(pp)
  pp
}
