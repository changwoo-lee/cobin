#' Random variate generation for micobin (mixture of continuous binomial) distribution
#'
#' Micobin distribution with natural parameter \eqn{\theta} and dispersion \eqn{psi}, denoted as \eqn{micobin(\theta, \psi)}, is defined as a dispersion mixture of cobin:
#' \deqn{
#'   Y \sim micobin(\theta, \psi) \iff Y | \lambda \sim cobin(\theta, \lambda^{-1}), (\lambda-1) \sim negbin(2, \psi) 
#' }
#'
#' @param n integer, number of samples
#' @param theta scalar or length n vector, natural parameter
#' @param psi scalar or length n vector, between 0 and 1, dispersion parameter
#' @param r (Default 2) This should be always 2 to maintain interpretaton of psi. It is kept for future experiment purposes. 
#'
#' @returns random samples from \eqn{micobin(\theta,\psi)}.
#' @export
#'
#' @examples
#' 
#' hist(rmicobin(1000, 2, 1/3), freq = FALSE)
#' xgrid = seq(0, 1, length = 500)
#' lines(xgrid, dmicobin(xgrid, 2, 1/3))
#' 
rmicobin <- function(n, theta, psi, r = 2){
  # ensure length of r is either 1 or n
  if(length(r) != 1){
    stop("length of r must be 1")
  }
  if(length(psi) != 1 & length(psi) != n){
    stop("length of psi must be either 1 or n")
  }
  # ensure psi is between 0 and 1
  if(any(psi < 0) || any(psi > 1)){
    stop("psi must be between 0 and 1")
  }
  # ensure length of theta is either 1 or n
  if(length(theta) != 1 & length(theta) != n){
    stop("length of theta must be either 1 or n")
  }

  if(length(theta)==1){
    lambda = rnbinom(n, r, psi) + 1
    return(rcobin(n, rep(theta, n), lambda))
  }else{
    lambda = rnbinom(n, r, psi) + 1
    return(rcobin(n, theta, lambda))
  }
}

#' Density function of micobin (mixture of continuous binomial) distribution
#' 
#' Micobin distribution with natural parameter \eqn{\theta} and dispersion \eqn{psi}, denoted as \eqn{micobin(\theta, \psi)}, is defined as a dispersion mixture of cobin:
#' \deqn{
#'   Y \sim micobin(\theta, \psi) \iff Y | \lambda \sim cobin(\theta, \lambda^{-1}), (\lambda-1) \sim negbin(2, \psi) 
#' }
#' so that micobin density is a weighted sum of cobin density with negative binomial weights.
#'
#' @param x num (length n), between 0 and 1, evaluation point
#' @param theta scalar or length n vector, natural parameter
#' @param psi scalar or length n vector, between 0 and 1, dispersion parameter
#' @param r (Default 2) This should be always 2 to maintain interpretaton of psi. It is kept for future experiment purposes.
#' @param log logical (Default FALSE), if TRUE, return log density
#' @param l_max integer (Default 70), upper bound of lambda.
#'
#' @returns density of \eqn{micobin(\theta, \psi)}
#' @export
#'
#' @examples
#' hist(rcobin(1000, 2, 3), freq = FALSE)
#' xgrid = seq(0, 1, length = 500)
#' lines(xgrid, dcobin(xgrid, 2, 3))
#' 
dmicobin <- function(x, theta, psi, r = 2, log = FALSE, l_max = 70){
  n = length(x)
#  if(length(psi) != 1){
#    stop("psi must be scalar")
#  }
  if(length(theta) == 1) theta = rep(theta, length(x))
  if(length(psi) == 1) psi = rep(psi, length(x))
  stopifnot("length of theta should be either 1 or length(x)" = (length(theta) == length(x)))
  stopifnot("length of psi should be either 1 or length(x)" = (length(psi) == length(x)))
  
  logdensity_summand = matrix(-Inf, l_max, n)
#  logweight_summand = matrix(-Inf, l_max, n)
 # for(l in 1:l_max) logweight_summand[l,] = dnbinom(l - 1, r, psi, log = T)
  logconst = pnbinom(l_max-1, r, psi, log.p = TRUE)
  maxerror = (1-pnbinom(l_max-1, r, psi, log.p = FALSE))
  if(any(maxerror > 0.0001)) warning("psi is small, so that deviation from truncated and untruncated may be large")
  for(l in 1:l_max){
    logdensity_summand[l,] = dnbinom(l - 1, r, psi, log = TRUE) - logconst + dcobin(x, theta, l, log = TRUE)
  }
  logdensity = matrixStats::colLogSumExps(logdensity_summand)
  if(log){
    return(logdensity)
  }else{
    return(exp(logdensity))
  }
}

#' Cumulative distribution function of micobin (mixture of continuous binomial) distribution
#' 
#' Micobin distribution with natural parameter \eqn{\theta} and dispersion \eqn{psi}, denoted as \eqn{micobin(\theta, \psi)}, is defined as a dispersion mixture of cobin:
#' \deqn{
#'   Y \sim micobin(\theta, \psi) \iff Y | \lambda \sim cobin(\theta, \lambda^{-1}), (\lambda-1) \sim negbin(2, \psi) 
#' }
#' so that micobin cdf is a weighted sum of cobin cdf with negative binomial weights.
#'
#' @param q num (length n), between 0 and 1, evaluation point
#' @param theta scalar, natural parameter
#' @param psi scalar, dispersion parameter
#' @param r (Default 2) This should be always 2 to maintain interpretaton of psi. It is kept for future experiment purposes.
#' @param l_max integer (Default 70), upper bound of lambda.
#'
#' @returns c.d.f. of \eqn{micobin(\theta, \psi)}
#' @export
#'
#' @examples
#' 
#' \donttest{
#' xgrid = seq(0, 1, length = 500)
#' out = pmicobin(xgrid, 1, 1/2)
#' plot(ecdf(rmicobin(10000, 1, 1/2)))
#' lines(xgrid, out, col = 2)
#'}
pmicobin <- function(q, theta, psi, r = 2, l_max = 70){
  if(length(psi) != 1){
    stop("psi must be scalar")
  }
  if(length(theta) != 1){
    stop("theta must be scalar")
  }
  cdf_summand = matrix(0, l_max, length(q))
  error = (1-pnbinom(l_max-1, r, psi, log.p = FALSE))
  if(error > 0.0001) warning("l_max may be too small")
  for(l in 1:l_max){
    cdf_summand[l,] = dnbinom(l - 1, r, psi)*pcobin(q, theta, l)
  }
  cdf = colSums(cdf_summand)
  return(cdf)
}

emicobin <- function(eta){
  bftprime(eta)
}

vmicobin <- function(eta, psi){
  bftprimeprime(eta)*psi
}

