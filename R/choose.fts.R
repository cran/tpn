
choose.fts<-function(y, criteria="AIC")
{
	if(any(y<=0)) stop("y must be positive")
	if(!any(criteria == c("AIC", "BIC"))) 
            stop("criteria is not recognized")
	aux.n =est.fts(y, dist="norm")
	aux.la=est.fts(y, dist="laplace")
	aux.c =est.fts(y, dist="cauchy")
	aux.lo=est.fts(y, dist="logis")
	if (grepl("Error", aux.n)[1] & grepl("Error", 
        aux.la)[1] & grepl("Error", aux.c)[1] & 
        grepl("Error", aux.lo)[1]) {
        stop("estimation problem in all the considered models")}
	index <- c()
    maxi = c()
    if (criteria == "AIC") {
        if (grepl("Error", aux.n)[1] == FALSE) {
            index <- c(index, 1)
            maxi <- c(maxi, aux.n$AIC)
        }
        if (grepl("Error", aux.la)[1] == FALSE) {
            index <- c(index, 2)
            maxi <- c(maxi, aux.la$AIC)
        }
        if (grepl("Error", aux.c)[1] == FALSE) {
            index <- c(index, 3)
            maxi <- c(maxi, aux.c$AIC)
        }
        if (grepl("Error", aux.lo)[1] == FALSE) {
            index <- c(index, 4)
            maxi <- c(maxi, aux.lo$AIC)
        }
    }
    if (criteria == "BIC") {
        if (grepl("Error", aux.n)[1] == FALSE) {
            index <- c(index, 1)
            maxi <- c(maxi, aux.n$BIC)
        }
        if (grepl("Error", aux.la)[1] == FALSE) {
            index <- c(index, 2)
            maxi <- c(maxi, aux.la$BIC)
        }
        if (grepl("Error", aux.c)[1] == FALSE) {
            index <- c(index, 3)
            maxi <- c(maxi, aux.c$BIC)
        }
        if (grepl("Error", aux.lo)[1] == FALSE) {
            index <- c(index, 4)
            maxi <- c(maxi, aux.lo$BIC)
        }
    }
	names(maxi)<-c("normal","laplace","cauchy","logistic")[index] 
	index <- index[which.min(maxi)]
	aux <- switch(index, "1" = aux.n, "2" = aux.la, 
        "3" = aux.c, "4" = aux.lo)
	selected <- aux$dist
	if(criteria=="AIC") res=list(AIC = maxi, selected = selected, estimate = aux$estimate,
		conv=aux$conv, logLik=aux$logLik, AIC=aux$AIC, BIC=aux$BIC)
	if(criteria=="BIC") res=list(BIC = maxi, selected = selected, estimate = aux$estimate,
		conv=aux$conv, logLik=aux$logLik, AIC=aux$AIC, BIC=aux$BIC)
	res
}
