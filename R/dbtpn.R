dbtpn <-
function(x, sigma, lambda, eta, log=FALSE)
{
    if (is.null(x)) 
        stop("x must be specified")
    if (is.null(sigma)) 
        stop("sigma must be specified")
    if (is.null(lambda)) 
        stop("lambda must be specified")
    if (is.null(eta)) 
        stop("eta must be specified")
    if (sigma <= 0) 
        stop("sigma must be positive")
lf=-log(2)-log(sigma)-pnorm(lambda, log.p=TRUE)+dnorm(x*sqrt(1+eta^2)/(sigma*(sqrt(1+eta^2)-eta))-lambda,log=TRUE)*ifelse(x>=0,1,0)+dnorm(-x*sqrt(1+eta^2)/(sigma*(sqrt(1+eta^2)+eta))-lambda,log=TRUE)*ifelse(x<0,1,0)
    if (!log) 
        lf = exp(lf)
    lf
}
