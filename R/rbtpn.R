rbtpn <-
function(n, sigma=1, lambda=0, eta=0)
{
    if (is.null(n)) 
        stop("sample size must be specified")
    if (is.null(sigma)) 
        stop("sigma must be specified")
    if (is.null(lambda)) 
        stop("lambda must be specified")
    if (is.null(eta)) 
        stop("eta must be specified")
    if (sigma <= 0) 
        stop("sigma must be positive")
    if (n <= 0 | round(n) != n) 
        stop("sample size must be a positive integer")
t=rtpn(n, sigma, lambda)
xi=eta/sqrt(1+eta^2)
z=ifelse(runif(n)<=(1+xi)/2, -1-xi, 1-xi)
t*z
}
