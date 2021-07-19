dtpn <-
function(x, sigma, lambda, log=FALSE)
{
if(is.null(x))
    stop("x must be specified")
  if(is.null(sigma))
    stop("sigma must be specified")
  if(is.null(lambda))
    stop("lambda must be specified")
  if(sigma<=0)
    stop("sigma must be positive")
  if(any(x<=0))
    stop("x's must be positive")
  lf=-log(sigma)-pnorm(lambda, log.p=TRUE)+dnorm(x/sigma-lambda, log=TRUE)
  if(!log) lf=exp(lf)
  lf
}
