\name{est.stpn}
\alias{est.stpn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter estimation for the stpn
}
\description{
Perform the parameter estimation for the slash truncated positive normal (stpn) discussed in Gomez, Gallardo and Santoro (2021)
based on the EM algorithm. Estimated errors are computed based on the Louis method to approximate the hessian matrix.
}
\usage{
est.stpn(y, sigma0=NULL, lambda0=NULL, q0=NULL, prec = 0.001, max.iter = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
the response vector. All the values must be positive.
}
  \item{sigma0, lambda0, q0}{
initial values for the EM algorithm for sigma, lambda and q. If they are omitted, 
by default sigma0 is defined as the root of the mean of the y^2, lambda as 0 and q as 3.
}
  \item{prec}{
the precision defined for each parameter. By default is 0.001.
}
  \item{max.iter}{
the maximum iterations for the EM algorithm. By default is 1000.
}
}
\details{
A variable have stpn distribution with parameters sigma>0, lambda in R and q>0 if its probability density
function can be written as

f(y; sigma, lambda, q) = int_0^1 (t^(1/q)*sigma)*phi(y*t^(1/q)*sigma-lambda)dt,

where phi() denotes the density function for the standard normal distribution.
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
Gomez, H., Gallardo, D.I., Santoro, K. (2021) Slash Truncation Positive Normal Distribution: with
application using the EM algorithm. Symmetry, 13, 2164.
}
\author{
Gallardo, D.I. and Gomez, H.J.
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
y=rstpn(n=100,sigma=10,lambda=1,q=2)
est.stpn(y)
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
