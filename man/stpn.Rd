\name{stpn}
\alias{dstpn}
\alias{pstpn}
\alias{rstpn}
\title{Slash truncated positive normal
}
\description{
Density, distribution function and random generation for the slash truncated positive normal (stpn) discussed in Gomez, Gallardo and Santoro (2021).
}
\usage{
dstpn(x, sigma, lambda, q, log = FALSE)
pstpn(x, sigma, lambda, q, lower.tail=TRUE, log=FALSE)
rstpn(n, sigma, lambda, q)
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
  \item{q}{
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
Random generation is based on the stochastic representation of the model, i.e., the quotient between a tpn
(see Gomez et al. 2018) and a beta random variable. 
}
\value{
dstpn gives the density, pstpn gives the distribution function and rstpn generates random deviates.

The length of the result is determined by n for rstpn, and is the maximum of the lengths of the numerical arguments for the other functions.

The numerical arguments other than n are recycled to the length of the result. Only the first elements of the logical arguments are used.

}
\references{
Gomez, H., Gallardo, D.I., Santoro, K. (2021) Slash Truncation Positive Normal Distribution: with
application using the EM algorithm. Submitted 

Gomez, H.J., Olmos, N.M., Varela, H., Bolfarine, H. (2018). Inference for a truncated positive normal
distribution. Applied Mathemetical Journal of Chinese Universities, 33, 163-176.
}
\author{
Gallardo, D.I. and Gomez, H.J.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%\code{\link{stpn}}
%}
\examples{
dstpn(c(1,2), sigma=1, lambda=-1, q=2)
pstpn(c(1,2), sigma=1, lambda=-1, q=2)
rstpn(n=10, sigma=1, lambda=-1, q=2)
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
