\name{est.btpn}
\alias{est.btpn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter estimation for the btpn
}
\description{
Perform the parameter estimation for the bimodal truncated positive normal (btpn) discussed in Gomez et al. (2021). 
Estimated errors are computed based on the hessian matrix.
}
\usage{
est.btpn(y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
the response vector. All the values must be positive.
}
}
\details{
A variable have btpn distribution with parameters sigma>0, lambda in R and eta in R if its probability density
function can be written as

f(y; sigma, lambda, q) = phi(x/(sigma*(1+epsilon))+lambda)/(2*sigma*Phi(lambda)), y<0
                       = phi(x/(sigma*(1-epsilon))-lambda)/(2*sigma*Phi(lambda)), y>=0
                         
where epsilon=eta/sqrt(1+eta^2) and phi() and Phi() denote the density function and the cumulative distribution 
function for the standard normal distribution, respectively.
}
\value{
A list with the following components
\item{estimate}{A matrix with the estimates and standard errors}
\item{iter}{Iterations in which the convergence were attached.}
\item{logLik}{log-likelihood function evaluated in the estimated parameters.}
\item{AIC}{Akaike's criterion.}
\item{BIC}{Schwartz's criterion.}
}
\references{
Gomez, H.J., Caimanque, W., Concha, M., Gomez, Y.M., Gallardo, D.I.  (2021) Bimodal Truncation Positive Normal 
Distribution. Submitted 
}
\author{
Gallardo, D.I., Gomez, H.J. and Gomez, Y.M.
}
\note{
A warning is presented if the estimated hessian matrix is not invertible.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%\code{\link{est.tpn}}
%}
\examples{
set.seed(2021)
y=rbtpn(n=100,sigma=10,lambda=1,eta=1.5)
est.btpn(y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.