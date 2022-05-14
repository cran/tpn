\name{ftp}
\alias{dftp}
\alias{pftp}
\alias{qftp}
\alias{rftp}
\title{Flexible truncated positive normal
}
\description{
Density, distribution function and random generation for the flexible truncated positive (ftp) class discussed in 
Gomez et al. (2022).
}
\usage{
dftp(x, sigma, lambda, dist="norm", log = FALSE)
pftp(x, sigma, lambda, dist="norm", lower.tail=TRUE, log.p=FALSE)
qftp(p, sigma, lambda, dist="norm")
rftp(n, sigma, lambda, dist="norm")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector of quantiles
}
  \item{p}{
vector of probabilities
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
  \item{dist}{
standard symmetrical distribution. Avaliable options: norm (default), logis, 
cauchy and laplace.
}
  \item{log, log.p}{
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
dftp gives the density, pftp gives the distribution function, qftp gives the quantile function and rftp generates random deviates.

The length of the result is determined by n for rbtpn, and is the maximum of the lengths of the numerical arguments for the other functions.

The numerical arguments other than n are recycled to the length of the result. Only the first elements of the logical arguments are used.

A variable have ftp distribution with parameters \eqn{\sigma>0} and \eqn{\lambda \in} R if its probability density
function can be written as
\deqn{
f(y; \sigma, \lambda, q) = \frac{g_0(\frac{y}{\sigma}-\lambda)}{\sigma G_0(\lambda)}, y>0,
}
where \eqn{g_0(\cdot)} and \eqn{G_0(\cdot)} denote the pdf and cdf for the specified distribution.
The case where \eqn{g_0(\cdot)} and \eqn{G_0(\cdot)} are from the standard normal model is known as the truncated positive normal model discussed
in Gomez et al. (2018).

}
\references{
Gomez, H.J., Gomez, H.W., Santoro, K.I., Venegas, O., Gallardo, D.I. (2022). A Family of Truncation Positive Distributions.
Submitted.

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
dftp(c(1,2), sigma=1, lambda=1, dist="logis")
pftp(c(1,2), sigma=1, lambda=1, dist="logis")
rftp(n=10, sigma=1, lambda=1, dist="logis")
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
