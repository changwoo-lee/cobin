% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/micobindist.R
\name{dmicobin}
\alias{dmicobin}
\title{Density function of micobin (mixture of continuous binomial) distribution}
\usage{
dmicobin(x, theta, psi, r = 2, log = FALSE, l_max = 70)
}
\arguments{
\item{x}{num (length n), between 0 and 1, evaluation point}

\item{theta}{scalar or length n vector, natural parameter}

\item{psi}{scalar or length n vector, between 0 and 1, dispersion parameter}

\item{r}{(Default 2) This should be always 2 to maintain interpretaton of psi. It is kept for future experiment purposes.}

\item{log}{logical (Default FALSE), if TRUE, return log density}

\item{l_max}{integer (Default 70), upper bound of lambda.}
}
\value{
density of \eqn{micobin(\theta, \psi)}
}
\description{
Micobin distribution with natural parameter \eqn{\theta} and dispersion \eqn{psi}, denoted as \eqn{micobin(\theta, \psi)}, is defined as a dispersion mixture of cobin:
\deqn{
  Y \sim micobin(\theta, \psi) \iff Y | \lambda \sim cobin(\theta, \lambda^{-1}), (\lambda-1) \sim negbin(2, \psi) 
}
so that micobin density is a weighted sum of cobin density with negative binomial weights.
}
\examples{
\dontrun{
hist(rcobin(1000, 2, 3), freq = FALSE)
xgrid = seq(0, 1, length = 500)
lines(xgrid, dcobin(xgrid, 2, 3))
}
}
