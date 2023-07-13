rutpn<-function(n, sigma=1, lambda=0, type=1)
{
if(type!=1 & type!=2) stop("type must be 1 or 2")
	if (is.null(n)) 
        stop("sample size must be specified")
    if (is.null(sigma)) 
        stop("sigma must be specified")
    if (is.null(lambda)) 
        stop("lambda must be specified")
	if(any(sigma<=0)) stop("sigma must be positive")
	if(type!=1 & type!=2) stop("type must be 1 or 2")
    if (n <= 0 | round(n) != n) 
        stop("sample size must be a positive integer")
	x<-rtpn(n, sigma, lambda)
	y<-1/(1+x)
	if(type==2) y<-1-y
	y 
}