putpn<-function (x, sigma=1, lambda=0, type=1, lower.tail = TRUE, log = FALSE) 
{
	if(type!=1 & type!=2 & type!=3 & type!=4) stop("type must be 1, 2, 3 or 4")
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
    if(type==1) lf = exp(-pnorm(lambda,log.p=TRUE)+log1p(-pnorm((1-x)/(sigma*x)-lambda)))
    if(type==2) lf = exp(-pnorm(lambda,log.p=TRUE)+log(-1+pnorm(lambda)+pnorm(x/(sigma*(1-x))-lambda)))
    if(type==3) lf = exp(-pnorm(lambda,log.p=TRUE)+pnorm(log(x)/sigma+lambda, log.p=TRUE))
    if(type==4) lf = exp(log1p(-exp(pnorm(log1p(-x)/sigma+lambda,log.p=TRUE)-pnorm(lambda,log.p=TRUE))))
    if (!lower.tail) 
        lf = 1 - lf
    if (log) 
        lf = log(lf)
    lf
}