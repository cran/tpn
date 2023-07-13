putpn<-function (x, sigma=1, lambda=0, type=1, lower.tail = TRUE, log = FALSE) 
{
	if(type!=1 & type!=2) stop("type must be 1 or 2")
    if (is.null(x)) 
        stop("x must be specified")
    if (is.null(sigma)) 
        stop("sigma must be specified")
    if (is.null(lambda)) 
        stop("lambda must be specified")
    if (sigma <= 0) 
        stop("sigma must be positive")
    if (any(x <= 0)) 
        stop("x's must be positive")
if(type!=1 & type!=2) stop("type must be 1 or 2")
    lf = exp(-pnorm(lambda,log.p=TRUE)+log1p(-pnorm((1-x)/(sigma*x)-lambda)))
	if(type==2) lf=exp(-pnorm(lambda,log.p=TRUE)+log(-1+pnorm(lambda)+pnorm(x/(sigma*(1-x))-lambda)))
    if (!lower.tail) 
        lf = 1 - lf
    if (log) 
        lf = log(lf)
    lf
}