rutpn<-function(n, sigma=1, lambda=0, type=1)
{
if(type!=1 & type!=2 & type!=3 & type!=4) stop("type must be 1, 2, 3 or 4")
	if (is.null(n)) 
        stop("sample size must be specified")
    if (is.null(sigma)) 
        stop("sigma must be specified")
    if (is.null(lambda)) 
        stop("lambda must be specified")
	if(any(sigma<=0)) stop("sigma must be positive")
    if (n <= 0 | round(n) != n) 
        stop("sample size must be a positive integer")
	x<-rtpn(n, sigma, lambda)
	if(type==1) y<-1/(1+x)
	if(type==2) y<-x/(1+x)
	if(type==3) y<-exp(-x)
	if(type==4) y<-1-exp(-x)
	y 
}