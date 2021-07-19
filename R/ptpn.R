ptpn <-
function(x, sigma, lambda, lower.tail=TRUE, log=FALSE)
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
  lf=1-exp(pnorm(x/sigma-lambda, log.p=TRUE, lower.tail = FALSE)-pnorm(lambda, log.p=TRUE))
  if(!lower.tail) lf=1-lf
  if(log) lf=log(lf)
  lf
}
