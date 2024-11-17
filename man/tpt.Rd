\name{tpt}
\alias{dtpt}
\alias{ptpt}
\alias{rtpt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Truncated positive t
}
\description{
Density, distribution function and random generation for the slash truncated positive student-t (tpt) distribution.
}
\usage{
dtpt(x, sigma, lambda, nu, log = FALSE)
ptpt(x, sigma, lambda, nu, lower.tail=TRUE, log=FALSE)
rtpt(n, sigma, lambda, nu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector of quantiles
}
  \item{n}{
number of observations
}
  \item{sigma}{
scale parameter for the distribution
}
  \item{lambda}{
shape parameter for the distribution
}
  \item{nu}{
nu parameter for the distribution
}
  \item{log}{
logical; if TRUE, probabilities p are given as log(p).
}
\item{lower.tail}{
logical; if TRUE (default), probabilities are P[X <= x] otherwise, P[X > x].
}
}
\details{
Random generation is based on the inverse transformation method. 
}
\value{
dtpt gives the density, ptpt gives the distribution function and rtpt generates random deviates.

The length of the result is determined by n for rtpt, and is the maximum of the lengths of the numerical arguments for the other functions.

The numerical arguments other than n are recycled to the length of the result. Only the first elements of the logical arguments are used.

A variable have tpt distribution with parameters \eqn{\sigma>0}, \eqn{\lambda \in} R and \eqn{\nu>0} if its probability density
function can be written as
\deqn{
f(y; \sigma, \lambda, q) = \frac{t_\nu\left(\frac{y}{\sigma}-\lambda\right)}{\sigma T_\nu(\lambda)}, y>0,
}
where \eqn{t_\nu(\cdot)} and \eqn{T_\nu(\cdot)} denote the density and cumultative distribution functions for the standard t distribution with \eqn{\nu} degrees of freedom.

}
%\references{
%Gomez, H.J., Olmos, N.M., Varela, H., Bolfarine, H. (2018). Inference for a truncated positive normal
%distribution. Applied Mathemetical Journal of Chinese Universities, 33, 163-176.
%}
\author{
Gallardo, D.I. and Gomez, H.J.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%\code{\link{tpn}}
%}
\examples{
dtpt(c(1,2), sigma=1, lambda=-1, nu=5)
ptpt(c(1,2), sigma=1, lambda=-1, nu=5)
rtpt(n=10, sigma=1, lambda=-1, nu=5)
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
