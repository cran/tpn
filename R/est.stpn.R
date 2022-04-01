est.stpn <-
function(y, sigma0=NULL, lambda0=NULL, q0=NULL, prec=1e-3,max.iter=1000)
{
if(is.null(y)) stop("y must be defined")
if(prec>0.1) stop("defined precision is too high")
if(max.iter<=0 | round(max.iter)!=max.iter) stop("iterations must be a positive integer")
if(is.null(sigma0)) sigma0=sqrt(mean(y^2))
if(is.null(lambda0)) lambda0=0
if(is.null(q0)) q0=3
if(sigma0<=0) stop("sigma0 must be positive")
if(q0<=0) stop("q0 must be positive")
EM.algorithm<-function(y,sigma.last,lambda.last,q.last,prec=1e-3,max.iter=1000){
C=function(y,theta,j=0,k=0){
aux.C=function(x,y,j,k,theta){
      sigma=theta[1]; lambda=theta[2]; q=theta[3];ll=x^(q+j)*log(x)^k*dnorm(x*y/sigma-lambda);ll}
integrate(aux.C,lower=0,upper=1,j=j,k=k,y=y,theta=theta)$value}
passo.E=function(y,theta){
C00<-sapply(X=y,FUN=C, theta=theta,j=0,k=0);C10<-sapply(X=y,FUN=C, theta=theta,j=1,k=0)
C20<-sapply(X=y,FUN=C, theta=theta,j=2,k=0);C01<-sapply(X=y,FUN=C, theta=theta,j=0,k=1)
u=C10/C00; u2=C20/C00; u10=C01/C00
return(list(u=u,u2=u2,u10=u10))}
passo.CM1=function(y,lambda,u,u2){
A=lambda*mean(y*u);B=mean(y^2*u2)
(-A+sqrt(A^2+4*B))/2}
equation.lambda<-function(lambda,sigma,y,u){
ll<--pnorm(lambda, log.p=TRUE)+y*u*lambda/sigma-lambda^2/2;-sum(ll)}
passo.CM2=function(y,lambda,sigma,u,u2){
optim(lambda, equation.lambda, sigma=sigma, y=y, u=u, method="BFGS")$par}
passo.CM3=function(y,u10){
-length(y)/(sum(u10))}
dif=10;i=1
while(dif>prec & i<=max.iter & q.last<40){
aux1=passo.E(y,c(sigma.last,lambda.last,q.last))
u=aux1$u;u2=aux1$u2;u10=aux1$u10
sigma.new=passo.CM1(y,lambda.last,u, u2)
lambda.new=passo.CM2(y,lambda.last,sigma.new, u)
q.new=passo.CM3(y,u10)
dif=max(abs(c(sigma.last-sigma.new,lambda.last-lambda.new,q.last-q.new)))
i=i+1;sigma.last=sigma.new;lambda.last=lambda.new;q.last=q.new}
par=c(sigma.last,lambda.last,q.last); names(par)=c("sigma", "lambda", "q")
conv=ifelse(i<max.iter, 0, 1) 
res<-list(par=par, iter=i-1, conv=conv)
if(sum(c(q.last>=20,i-1==max.iter))==1)
{if(q.last>=20) res$warnings="q is too large. It is suggested the tpn instead of the stpn model"
if(i-1==max.iter) res$warnings="the iterations limit had been reached without convergence"}
if(sum(c(q.last>=10,i-1==max.iter))==2)
{
res$warnings="q is too large. It is suggested the tpn instead of the stpn model"
res$warnings[2]="the iterations limit had been reached"
}
res
}
llike.STPN<-function(theta, y)
{
sigma=theta[1];lambda=theta[2];q=theta[3]
aux.int=function(x, y, lambda, sigma, q){
exp(q*log(x)+dnorm(x*y/sigma-lambda,log=TRUE))}
int=c()
for(i in 1:length(y)){
int[i]=integrate(aux.int, lower=0, upper=1, y=y[i], lambda=lambda, sigma=sigma, q=q)$value}
-sum(log(q)-log(sigma)-pnorm(lambda, log.p=TRUE)+log(int))
}
Louis.STPN=function(theta,y)
{sigma=theta[1];lambda=theta[2];q=theta[3]
C=function(y,theta,j=0,k=0){
aux.C=function(x,y,j,k,theta){
      sigma=theta[1];lambda=theta[2];q=theta[3]
      ll=x^(q+j)*log(x)^k*dnorm(x*y/sigma-lambda);ll}
integrate(aux.C,lower=0,upper=1,j=j,k=k,y=y,theta=theta)$value}
C00<-sapply(X=y,FUN=C, theta=theta,j=0,k=0);C10<-sapply(X=y,FUN=C, theta=theta,j=1,k=0)
C20<-sapply(X=y,FUN=C, theta=theta,j=2,k=0);C30<-sapply(X=y,FUN=C, theta=theta,j=3,k=0)
C40<-sapply(X=y,FUN=C, theta=theta,j=4,k=0);C01<-sapply(X=y,FUN=C, theta=theta,j=0,k=1)
C02<-sapply(X=y,FUN=C, theta=theta,j=0,k=2);C11<-sapply(X=y,FUN=C, theta=theta,j=1,k=1)
C21<-sapply(X=y,FUN=C, theta=theta,j=2,k=1);u=C10/C00;u2=C20/C00;u10=C01/C00;u3=C30/C00
u4=C40/C00;u11=C11/C00;u21=C21/C00;u20=C02/C00
B=matrix(0,ncol=3,nrow=3)
for(i in (1:length(y))){
z=exp(dnorm(lambda,log=TRUE)-pnorm(lambda,log.p=TRUE)) 
B[1,1]=B[1,1]+1/sigma^2-(3*y[i]^2*u2[i])/(sigma^4)+(2*y[i]*u[i]*lambda)/(sigma^3)
B[1,2]=B[2,1]=B[1,2]-(y[i]*u[i])/sigma^2
B[1,3]=B[3,1]=B[1,3]+0
B[2,2]=B[2,2]+z*(lambda+z)-1
B[2,3]=B[3,2]=B[2,3]+0
B[3,3]=B[3,3]-1/q^2}
D=matrix(0,ncol=3,nrow=3)
for(i in (1:length(y))){
z=exp(dnorm(lambda,log=TRUE)-pnorm(lambda,log.p=TRUE)) 
D[1,1]=D[1,1]+1/sigma^2-(2*y[i]^2*u2[i])/(sigma^4)+(2*y[i]*u[i]*lambda)/(sigma^3)+(y[i]^4*u4[i])/(sigma^6)-(2*y[i]^3*u3[i]*lambda)/(sigma^5)+(y[i]^2*u2[i]*lambda^2)/(sigma^4)
D[1,2]=D[2,1]=D[1,2]+z/sigma+((y[i]*u[i])/(sigma^2))*(lambda^2-1+lambda*z)-((y[i]^2*u2[i])/(sigma^3))*(z+lambda)+(y[i]^3*u3[i])/(sigma^4)-(y[i]^2*u2[i]*lambda)/(sigma^3)+lambda/sigma
D[1,3]=D[3,1]=D[1,3]-1/(sigma*q)-(u10[i])/(sigma)+(y[i]^2*u2[i])/(q*sigma^3)+(y[i]^2*u21[i])/(sigma^3)-(y[i]*u[i]*lambda)/(q*sigma^2)-(lambda*y[i]*u11[i])/(sigma^2)
D[2,2]=D[2,2]+z^2+(y[i]^2*u2[i])/(sigma^2)+lambda^2-(2*z*y[i]*u[i])/(sigma)+2*z*lambda-(2*y[i]*u[i]*lambda)/(sigma)
D[2,3]=D[3,2]=D[2,3]-z/q-z*u10[i]+(y[i]*u[i])/(q*sigma)+(y[i]*u11[i])/(sigma)-lambda/q-lambda*u10[i]
D[3,3]=D[3,3]+1/q^2+(2*u10[i])/(q)+u20[i]}
F=matrix(0,ncol=3,nrow=3)
for(i in 1:length(y)){
for(j in 1:length(y)){
if(i!=j)
{z=exp(dnorm(lambda,log=TRUE)-pnorm(lambda,log.p=TRUE)) 
F[1,1]=F[1,1]+1/sigma^2+(y[i]^2*u2[i]*y[j]^2*u2[j])/(sigma^6)+(y[i]*u[i]*y[j]*u[j]*lambda^2)/(sigma^4)-(y[i]^2*u2[i])/(sigma^4)-(y[j]^2*u2[j])/(sigma^4)-(lambda*y[j]^2*u2[j]*y[i]*u[i])/(sigma^5)-(lambda*y[i]^2*u2[i]*y[j]*u[j])/(sigma^5)+(y[i]*u[i]*lambda)/(sigma^3)+(y[j]*u[j]*lambda)/(sigma^3)                                       
F[1,2]=F[1,2]+z/sigma-(y[j]*u[j])/(sigma^2)+lambda/sigma-(y[i]^2*u2[i]*z)/(sigma^3)+(y[i]^2*y[j]*u2[i]*u[j])/(sigma^4)-(lambda*y[i]^2*u2[i])/(sigma^3)+(y[i]*u[i]*z*lambda)/(sigma^2)-(lambda*y[i]*u[i]*y[j]*u[j])/(sigma^3)+(lambda^2*y[i]*u[i])/(sigma^2)
F[2,1]=F[2,1]+z/sigma-(y[i]*u[i])/(sigma^2)+lambda/sigma-(y[j]^2*u2[j]*z)/(sigma^3)+(y[j]^2*y[i]*u2[j]*u[i])/(sigma^4)-(lambda*y[j]^2*u2[j])/(sigma^3)+(y[j]*u[j]*z*lambda)/(sigma^2)-(lambda*y[j]*u[j]*y[i]*u[i])/(sigma^3)+(lambda^2*y[j]*u[j])/(sigma^2)
F[1,3]=F[1,3]-1/(sigma*q)+(y[i]^2*u2[i])/(q*sigma^3)-(y[i]*u[i]*lambda)/(q*sigma^2)-(u10[j])/sigma+(y[i]^2*u2[i]*u10[j])/(sigma^3)-(y[i]*u[i]*u10[j]*lambda)/(sigma^2)
F[3,1]=F[3,1]-1/(sigma*q)+(y[j]^2*u2[j])/(q*sigma^3)-(y[j]*u[j]*lambda)/(q*sigma^2)-(u10[i])/sigma+(y[j]^2*u2[j]*u10[i])/(sigma^3)-(y[j]*u[j]*u10[i]*lambda)/(sigma^2)
F[2,3]=F[2,3]-z/q-z*u10[j]+(y[i]*u[i])/(q*sigma)+(y[i]*u[i]*u10[j])/(sigma)-lambda/q-lambda*u10[j]
F[3,2]=F[3,2]-z/q-z*u10[i]+(y[j]*u[j])/(q*sigma)+(y[j]*u[j]*u10[i])/(sigma)-lambda/q-lambda*u10[i]
F[2,2]=F[2,2]+z^2+(y[i]*u[i]*y[j]*u[j])/(sigma^2)+lambda^2-(z*y[i]*u[i])/(sigma)-(z*y[j]*u[j])/(sigma)-(lambda*y[i]*u[i])/(sigma)-(lambda*y[j]*u[j])/(sigma)+2*lambda*z
F[3,3]=F[3,3]+1/q^2+u10[i]*u10[j]+u10[i]/q+u10[j]/q}}};B+D+F
}
aux=try(EM.algorithm(y,sigma0,lambda0,q0,prec=prec,max.iter=max.iter),silent=TRUE)
aux2=try(solve(-Louis.STPN(aux$par,y)), silent=TRUE)
res<-c();res$warnings<-"Estimates can't be computed. Try with different initial values."
if(!grepl("Error",aux)[1])
{param=cbind(aux$par); colnames(param)<-c("estimate");rownames(param)<-c("sigma", "lambda", "q");llike=-llike.STPN(param[,1], y) 
if(!grepl("Error",aux2)[1])
{if(min(diag(aux2))>0)
{param=cbind(param, sqrt(diag(aux2)))
colnames(param)<-c("estimate", "s.e.")}}
res<-list(estimate=param, iter=aux$iter, conv=aux$conv, logLik=llike,
AIC=-2*llike+2*3, BIC=-2*llike+3*log(length(y)))
if(is.null(aux$warnings))
{if(ncol(param)==1) res$warnings="Standard errors can't be estimated"}
if(!is.null(aux$warnings))
{res$warnings=aux$warnings;
if(ncol(param)==1) res$warnings[length(res$warnings)+1]="Standard errors can't be estimated"}
}
res
}
