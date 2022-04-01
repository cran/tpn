\name{est.ftpn}
\alias{est.ftpn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter estimation for the ftpn
}
\description{
Perform the parameter estimation for the Flexible truncated positive normal (ftpn) discussed in Gomez et al. (2022)
based on maximum likelihood estimation. Estimated errors are computed based on the hessian matrix.
}
\usage{
est.ftpn(y, dist="norm")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
the response vector. All the values must be positive.
}
  \item{dist}{
standard symmetrical distribution. Avaliable options: norm (default), logis, 
cauchy and laplace.
}
}
\details{
A variable have ftpn distribution with parameters sigma>0 and lambda in R, if its probability density
function can be written as

f(y; sigma, lambda, q) = (sigma*G0(lambda))^(-1)*g0(y/sigma-lambda),

where g0() and G0() denote the density and cumultative distribution functions for the specified distribution.
The case where g0() and G0() are from the normal model is known as the truncated positive normal model discussed
in Gomez et al. (2018).
}
\value{
A list with the following components
\item{estimate}{A matrix with the estimates and standard errors}
\item{dist}{distribution specified}
\item{conv}{the code related to the convergence for the optim function. 0 if the convergence was attached.}
\item{logLik}{log-likelihood function evaluated in the estimated parameters.}
\item{AIC}{Akaike's criterion.}
\item{BIC}{Schwartz's criterion.}
}
\references{
Gomez, H.J., Gomez, H.W., Santoro, K.I., Venegas, O., Gallardo, D.I. (2022). A Family of Truncation Positive Distributions.
Submitted.

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
y=rftpn(n=100,sigma=10,lambda=1,dist="logis")
est.ftpn(y,dist="logis")
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
