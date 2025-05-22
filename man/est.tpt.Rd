\name{est.tpt}
\alias{est.tpt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter estimation for the tpt distribution
}
\description{
Perform the parameter estimation for the truncated positive Student's-t (tpt) distribution
based on maximum likelihood estimation. Estimated errors are computed based on the hessian matrix.
}
\usage{
est.tpt(y, x = NULL, q = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
the response vector. All the values must be positive.
}
\item{x}{
the covariates vector.
}
\item{q}{
quantile of the distribution to be modelled.
}
}
\details{
A variable have tpt distribution with parameters \eqn{\sigma>0}, \eqn{\lambda \in} R and \eqn{\nu>0} if its probability density
function can be written as
\deqn{
f(y; \sigma, \lambda, q) = \frac{t_\nu\left(\frac{y}{\sigma}-\lambda\right)}{\sigma T_\nu(\lambda)}, y>0,
}
where \eqn{t_\nu(\cdot)} and \eqn{T_\nu(\cdot)} denote the density and cumulative distribution functions for the standard t distribution with 
\eqn{\nu} degrees of freedom.
}
\value{
A list with the following components
\item{estimate}{A matrix with the estimates and standard errors}
\item{logLik}{log-likelihood function evaluated in the estimated parameters.}
\item{AIC}{Akaike's criterion.}
\item{BIC}{Schwartz's criterion.}
}
%\references{
%Gomez, H.J., Olmos, N.M., Varela, H., Bolfarine, H. (2018). Inference for a truncated positive normal
%distribution. Applied Mathemetical Journal of Chinese Universities, 33, 163-176.
%}
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
y=rtpt(n=100,sigma=10,lambda=1, nu=5)
est.tpt(y)
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
