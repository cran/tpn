\name{btpn}
\alias{dbtpn}
\alias{pbtpn}
\alias{rbtpn}
\title{Bimodal truncated positive normal
}
\description{
Density, distribution function and random generation for the bimodal truncated positive normal (btpn) discussed in 
Gomez et al. (2022).
}
\usage{
dbtpn(x, sigma, lambda, eta, log = FALSE)
pbtpn(x, sigma, lambda, eta, lower.tail=TRUE, log=FALSE)
rbtpn(n, sigma, lambda, eta)
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
  \item{eta}{
shape parameter for the distribution
}
  \item{log}{
logical; if TRUE, probabilities p are given as log(p).
}
\item{lower.tail}{
logical; if TRUE (default), probabilities are P[X <= x] otherwise, P[X > x].
}
}
\details{
Random generation is based on the stochastic representation of the model, i.e., the product between a tpn
(see Gomez et al. 2018) and a dichotomous variable assuming values \eqn{-(1+\epsilon)} and \eqn{1-\epsilon} with probabilities
\eqn{(1+\epsilon)/2} and \eqn{(1-\epsilon)/2}, respectively. 
}
\value{
dbtpn gives the density, pbtpn gives the distribution function and rbtpn generates random deviates.

The length of the result is determined by n for rbtpn, and is the maximum of the lengths of the numerical arguments for the other functions.

The numerical arguments other than n are recycled to the length of the result. Only the first elements of the logical arguments are used.

A variable have btpn distribution with parameters \eqn{\sigma>0, \lambda \in} R and \eqn{\eta \in} R if its probability density
function can be written as
\deqn{
f(y; \sigma, \lambda, q) = \frac{\phi\left(\frac{x}{\sigma(1+\epsilon)}+\lambda\right)}{2\sigma\Phi(\lambda)}, y<0,
}                         
and
\deqn{
f(y; \sigma, \lambda, q) = \frac{\phi\left(\frac{x}{\sigma(1-\epsilon)}-\lambda\right)}{2\sigma\Phi(\lambda)}, y\geq 0,
}
where \eqn{\epsilon=\eta/\sqrt{1+\eta^2}} and \eqn{\phi(\cdot)} and \eqn{\Phi(\cdot)} denote the probability density function and the cumulative distribution 
function for the standard normal distribution, respectively.
}
\references{
Gomez, H.J., Caimanque, W., Gomez, Y.M., Magalhaes, T.M., Concha, M., Gallardo, D.I. (2022) Bimodal Truncation Positive Normal 
Distribution. Symmetry, 14, 665.

Gomez, H.J., Olmos, N.M., Varela, H., Bolfarine, H. (2018). Inference for a truncated positive normal
distribution. Applied Mathemetical Journal of Chinese Universities, 33, 163-176.
}
\author{
Gallardo, D.I., Gomez, H.J. and Gomez, Y.M.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%\code{\link{stpn}}
%}
\examples{
dbtpn(c(1,2), sigma=1, lambda=-1, eta=2)
pbtpn(c(1,2), sigma=1, lambda=-1, eta=2)
rbtpn(n=10, sigma=1, lambda=-1, eta=2)
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
