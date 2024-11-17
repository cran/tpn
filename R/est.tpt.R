est.tpt <-
function(y)
{
llike.TPT<-function(theta, y, trans.param=FALSE)
{
sigma=theta[1]; lambda=theta[2]; nu=theta[3]; if(trans.param) {sigma=exp(theta[1]);nu=exp(theta[3])}
ll=-log(sigma)-pt(lambda, df=nu, log.p=TRUE)+dt(y/sigma-lambda, df=nu, log=TRUE);-sum(ll)
}
aux=optim(c(0,0,0), llike.TPT, y=y, trans.param=TRUE, method="BFGS", control=list(maxit=10000))
param=cbind(c(exp(aux$par[1]), aux$par[2], exp(aux$par[3])))
colnames(param)=c("estimate");se=try(solve2(hessian(llike.TPT, x0=param, y=y)),silent=TRUE)
llike=-aux$value
if(!grepl("Error",se)[1])
{if(min(diag(se))>0)
{param=cbind(param, sqrt(diag(se)))
colnames(param)<-c("estimate", "s.e.")}}
rownames(param)<-c("sigma", "lambda", "nu")
res=list(estimate=param, conv=aux$conv, logLik=llike, AIC=-2*llike+2*3, BIC=-2*llike+3*log(length(y)))
if(ncol(param)==1) res$warnings="Standard errors can't be estimated"
res
}
