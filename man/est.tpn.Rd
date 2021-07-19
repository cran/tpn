\name{est.tpn}
\alias{est.tpn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter estimation for the tpn
}
\description{
Perform the parameter estimation for the truncated positive normal (tpn) discussed in Gomez et al. (2018)
based on maximum likelihood estimation. Estimated errors are computed based on the hessian matrix.
}
\usage{
est.tpn(y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
the response vector. All the values must be positive.
}
}
\details{
A variable have tpn distribution with parameters sigma>0 and lambda in R, if its probability density
function can be written as

f(y; sigma, lambda, q) = (sigma*Phi(lambda))^(-1)*phi(y/sigma-lambda)dt,

where phi() and Phi() denote the density and cumultative distribution functions for the standard normal distribution.
}
\value{
A list with the following components
\item{estimate}{A matrix with the estimates and standard errors}
\item{logLik}{log-likelihood function evaluated in the estimated parameters.}
\item{AIC}{Akaike's criterion.}
\item{BIC}{Schwartz's criterion.}
}
\references{
Gomez, H.J., Olmos, N.M., Varela, H., Bolfarine, H. (2018). Inference for a truncated positive normal
distribution. Applied Mathemetical Journal of Chinese Universities, 33, 163-176.
}
\author{
Gallardo, D.I. and Gomez, H.J.
}
\note{
A warning is presented if the estimated hessian matrix is not invertible.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%\code{\link{est.stpn}}
%}
\examples{
set.seed(2021)
y=rtpn(n=100,sigma=10,lambda=1)
est.tpn(y)
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
