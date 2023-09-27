\name{est.utpn}
\alias{est.utpn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter estimation for the utpn model
}
\description{
Perform the parameter estimation for the unit truncated positive normal (utpn) type 1, 2, 3 or 4,
parameterized in terms of the quantile based on maximum likelihood estimation. 
Estimated errors are computed based on the hessian matrix.
}
\usage{
est.utpn(y, x=NULL, type=1, link="logit", q=0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
the response vector. All the values must be positive.
}
\item{x}{
the covariates vector.
}
\item{type}{
to distinguish the type of the utpn model: 1 (default), 2, 3 or 4.
}
\item{link}{
link function to be used for the covariates: logit (default).
}
\item{q}{
quantile of the distribution to be modelled.
}
}
%\details{
%A variable have tpn distribution with parameters \eqn{\sigma>0} and \eqn{\lambda \in} R if its probability density
%function can be written as
%\deqn{
%f(y; \sigma, \lambda, q) = \frac{\phi\left(\frac{y}{\sigma}-\lambda\right)}{\sigma \Phi(\lambda)}, y>0,
%}
%where \eqn{\phi(\cdot)} and \eqn{\Phi(\cdot)} denote the density and cumultative distribution functions for the standard normal distribution.
%}
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
Gallardo, D.I.
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
y=rutpn(n=100,sigma=10,lambda=1)
est.utpn(y)
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
