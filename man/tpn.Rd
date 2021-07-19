\name{tpn}
\alias{dtpn}
\alias{ptpn}
\alias{rtpn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Truncated positive normal
}
\description{
Density, distribution function and random generation for the slash truncated positive normal (stpn) discussed in Gomez, Gallardo and Santoro (2021).
}
\usage{
dtpn(x, sigma, lambda, log = FALSE)
ptpn(x, sigma, lambda, lower.tail=TRUE, log=FALSE)
rtpn(n, sigma, lambda)
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
dtpn gives the density, ptpn gives the distribution function and rtpn generates random deviates.

The length of the result is determined by n for rtpn, and is the maximum of the lengths of the numerical arguments for the other functions.

The numerical arguments other than n are recycled to the length of the result. Only the first elements of the logical arguments are used.

}
\references{
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
%\code{\link{tpn}}
%}
\examples{
dtpn(c(1,2), sigma=1, lambda=-1)
ptpn(c(1,2), sigma=1, lambda=-1)
rtpn(n=10, sigma=1, lambda=-1)
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
