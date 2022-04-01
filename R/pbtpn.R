pbtpn <-
function(x, sigma=1, lambda=0, eta=0, lower.tail = TRUE, log = FALSE)
{
    if (is.null(x)) 
        stop("x must be specified")
    if (is.null(sigma)) 
        stop("sigma must be specified")
    if (is.null(lambda)) 
        stop("lambda must be specified")
    if (sigma <= 0) 
        stop("sigma must be positive")
    epsilon = eta/(sqrt(1+eta^2))
    lf1 = 0.5*(1+epsilon)*(exp(-pnorm(lambda, log.p=TRUE))-exp(pnorm(-x/(sigma*(1+epsilon))-lambda, log.p=TRUE)-pnorm(lambda, log.p=TRUE)))
    lf2 = 1-0.5*(1-epsilon)*(exp(-pnorm(lambda, log.p=TRUE))-exp(pnorm(x/(sigma*(1-epsilon))-lambda, log.p=TRUE)-pnorm(lambda, log.p=TRUE)))
    lf3 = ifelse(x<=0,1,0)
    lf = lf1*lf3+lf2*(1-lf3)
    if (!lower.tail) 
        lf = 1 - lf
    if (log) 
        lf = log(lf)
    lf
}
