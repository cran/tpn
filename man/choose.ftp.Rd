\name{choose.ftp}
\alias{choose.ftp}
\title{
Choose a distribution in the flexible truncated positive class of models
}
\description{
Provide model selection for a given data set in the flexible truncated positive class of models
}
\usage{
choose.ftp(y, criteria = "AIC")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
positive vector of responses
}
  \item{criteria}{
model criteria for the selection: AIC (default) or BIC.
}
}
\details{
The function fits the truncated positive normal, truncated positive laplace, 
truncated positive Cauchy and truncated positive logistic models and select
the model which provides the lower criteria (AIC or BIC).
}
\value{
A list with the following components
\item{AIC}{a vector with the AIC for the different truncated positive fitted models:
	normal, laplace, cauchy and logistic.}
\item{selected}{the selected model}
\item{estimate}{the estimated for sigma and lambda and the respective standard errors (s.e.)}
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
Gallardo, D.I., Gomez, H.J. and Gomez, Y.M.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
set.seed(2021)
y=rftp(n=100,sigma=10,lambda=1,dist="logis")
choose.ftp(y)
}


