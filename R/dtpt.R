dtpt <-
function(x, sigma, lambda, nu, log=FALSE)
{
if(is.null(x))
    stop("x must be specified")
  if(is.null(sigma))
    stop("sigma must be specified")
  if(is.null(lambda))
    stop("lambda must be specified")
  if(is.null(nu))
    stop("nu must be specified")
  if(sigma<=0)
    stop("sigma must be positive")
  if(nu<=0)
    stop("nu must be positive")
  if(any(x<=0))
    stop("x's must be positive")
  lf=-log(sigma)-pt(lambda, df=nu, log.p=TRUE)+dt(x/sigma-lambda, df=nu, log=TRUE)
  if(!log) lf=exp(lf)
  lf
}
