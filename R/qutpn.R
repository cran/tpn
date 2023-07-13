qutpn<-function(p, sigma=1, lambda=0, type=1)
{
	if(type!=1 & type!=2) stop("type must be 1 or 2")
	 if (is.null(p)) 
        stop("p must be specified")
    if (is.null(sigma)) 
        stop("sigma must be specified")
    if (is.null(lambda)) 
        stop("lambda must be specified")
    if (sigma <= 0) 
        stop("sigma must be positive")
    if (any(p <= 0 | p>=1)) 
        stop("x's must be between 0 and 1")
	if(type!=1 & type!=2) stop("type must be 1 or 2")
	qq=(1+sigma*(lambda+qnorm(1-p*pnorm(lambda))))^(-1)
	if(type==2) qq=sigma*(qnorm(pnorm(lambda)*(p-1)+1))/(1+sigma*qnorm(pnorm(lambda)*(p-1)+1)+lambda)
	qq
}