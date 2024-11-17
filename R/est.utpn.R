est.utpn=function (y, x = NULL, type = 1, link = "logit", q = 0.5) 
{
    if (any(y <= 0 | y >= 1)) 
        stop("y must be between 0 and 1")
    if (!any(link == c("logit"))) 
        stop("link is not recognized")
    if (type != 1 & type != 2 & type != 3 & type != 4) 
        stop("type must be 1, 2, 3 or 4")
    if (is.null(x)) {
        x <- matrix(1, ncol = 1, nrow = length(y))
        colnames(x) <- "Intercept"
    }
    if (!is.matrix(x)) 
        stop("x must be a matrix")
    if (nrow(x) != length(y)) 
        stop("rows of x don't agree with the length of y")
    r = ncol(x)
    if (is.null(colnames(x))) 
        colnames(x) <- paste("var", 1:r, sep = "")
    if (r > 1) {
        llike.utpn <- function(theta, y, x, type = 1, link = "logit", 
            q = 0.5, t.param = FALSE) {
            beta = theta[1:r]
            lambda = theta[r + 1]
            dist <- switch(link, logit = "logis", probit = "norm", 
			loglog = "gumbel", cloglog = "gumbel2")
            g <- get(paste("d", dist, sep = ""), mode = "function")
            G <- get(paste("p", dist, sep = ""), mode = "function")
            rho = G(x %*% beta)
            if (type == 1) {
                log.sigma <- log1p(-rho) - log(rho) - log(lambda + 
                  qnorm(1 - q * pnorm(lambda)))
                ll <- -log.sigma - 2 * log(y) - pnorm(lambda, 
                  log.p = TRUE) + dnorm((1 - y) * exp(-log.sigma)/y - 
                  lambda, log = TRUE)
            }
            if (type == 2) {
                log.sigma <- log(rho) - log1p(-rho) - log(lambda + 
                  qnorm(1 - (1 - q) * pnorm(lambda)))
                ll <- -log.sigma - 2 * log1p(-y) - pnorm(lambda, 
                  log.p = TRUE) + dnorm(y * exp(-log.sigma)/(1 - 
                  y) - lambda, log = TRUE)
            }
            if (type == 3) {
                log.sigma <- log(-log(rho)) - log(lambda - qnorm(q * 
                  pnorm(lambda)))
                ll <- -log.sigma - log(y) - pnorm(lambda, log.p = TRUE) + 
                  dnorm(log(y) * exp(-log.sigma) + lambda, log = TRUE)
            }
            if (type == 4) {
                log.sigma <- log(-log1p(-rho)) - log(lambda - 
                  qnorm((1 - q) * pnorm(lambda)))
                ll <- -log.sigma - log1p(-y) - pnorm(lambda, 
                  log.p = TRUE) + dnorm(log1p(-y) * exp(-log.sigma) + 
                  lambda, log = TRUE)
            }
            -sum(ll)
        }
        type.aux = type
        aux <- optim(rep(0, r + 1), llike.utpn, y = y, x = x, 
            q = q, type = type.aux, link = link, method = "BFGS", 
            control = list(maxit = 10000))
        param = matrix(aux$par, ncol = 1)
        colnames(param) = c("estimate")
        llike = -aux$value
	se = try(solve(aux$hessian),silent=TRUE)
        if (!grepl("Error", se)[1]) {
            if (min(diag(se)) > 0) {
                param = cbind(param, sqrt(diag(se)))
                colnames(param) <- c("estimate", "s.e.")
            }
        }
        rownames(param) <- c(paste("beta", colnames(x), sep = ""), 
            "lambda")
        res = list(estimate = param, dist = paste("utpn", type, 
            sep = ""), quantile = q, conv = aux$conv, logLik = llike, 
            AIC = -2 * llike + 2 * nrow(param), BIC = -2 * llike + 
                nrow(param) * log(length(y)))
        if (ncol(param) == 1) 
            res$warnings = "Standard errors can't be estimated"
    }
    if (r == 1) {
        llike.utpn0 <- function(theta, y, type = 1, link = "logit", 
            q = 0.5, t.param = FALSE) {
            sigma = theta[1]
            lambda = theta[2]
            if (t.param) 
                sigma = exp(theta[1])
            dist <- switch(link, logit = "logis")
            g <- get(paste("d", dist, sep = ""), mode = "function")
            G <- get(paste("p", dist, sep = ""), mode = "function")
            if (type == 1) {
                log.sigma <- log(sigma)
                ll <- -log.sigma - 2 * log(y) - pnorm(lambda, 
                  log.p = TRUE) + dnorm((1 - y) * exp(-log.sigma)/y - 
                  lambda, log = TRUE)
            }
            if (type == 2) {
                log.sigma <- log(sigma)
                ll <- -log.sigma - 2 * log1p(-y) - pnorm(lambda, 
                  log.p = TRUE) + dnorm(y * exp(-log.sigma)/(1 - 
                  y) - lambda, log = TRUE)
            }
            if (type == 3) {
                log.sigma <- log(sigma)
                ll <- -log.sigma - log(y) - pnorm(lambda, log.p = TRUE) + 
                  dnorm(log(y) * exp(-log.sigma) + lambda, log = TRUE)
            }
            if (type == 4) {
                log.sigma <- log(sigma)
                ll <- -log.sigma - log1p(-y) - pnorm(lambda, 
                  log.p = TRUE) + dnorm(log1p(-y) * exp(-log.sigma) + 
                  lambda, log = TRUE)
            }
            -sum(ll)
        }
        type.aux = type
        aux <- optim(rep(0, 2), llike.utpn0, y = y, q = q, type = type.aux, 
            link = link, method = "BFGS", control = list(maxit = 10000), 
            t.param = TRUE)
        param = matrix(c(exp(aux$par[1]), aux$par[2]), ncol = 1)
        colnames(param) = c("estimate")
        llike = -aux$value
        se = try(skewMLRM::solve2(pracma::hessian(llike.utpn0, 
            x0 = param, y = y, type = type, link = link)), silent = TRUE)
        if (!grepl("Error", se)[1]) {
            if (min(diag(se)) > 0) {
                param = cbind(param, sqrt(diag(se)))
                colnames(param) <- c("estimate", "s.e.")
            }
        }
        rownames(param) <- c("sigma", "lambda")
        res = list(estimate = param, dist = paste("utpn", type, 
            sep = ""), conv = aux$conv, logLik = llike, AIC = -2 * 
            llike + 2 * nrow(param), BIC = -2 * llike + nrow(param) * 
            log(length(y)))
        if (ncol(param) == 1) 
            res$warnings = "Standard errors can't be estimated"
    }
    res
}