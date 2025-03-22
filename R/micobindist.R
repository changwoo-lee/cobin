#' Title
#'
#' @param n 
#' @param theta 
#' @param psi 
#' @param r 
#'
#' @returns
#' @export
#'
#' @examples
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

#' Title
#'
#' @param x 
#' @param theta 
#' @param psi 
#' @param r 
#' @param log 
#' @param l_max 
#'
#' @returns
#' @export
#'
#' @examples
dmicobin <- function(x, theta, psi, r = 2, log = FALSE, l_max = 70){
  if(length(psi) != 1){
    stop("psi must be scalar")
  }
  if(length(theta) == 1){
    theta = rep(theta, length(x))
  }
  # ensure legnth of theta is also n; if scalar, convert it to rep(theta, n)
  if(length(theta) != length(x)){
    stop("theta must be scalar or same as length of x")
  }
  logdensity_summand = matrix(-Inf, l_max, length(x))
  error = (1-pnbinom(l_max-1, r, psi, log = FALSE))
  if(error > 0.0001) warning("m_max is too small and density may be off by 0.0001. Increase l_max.")
  for(l in 1:l_max){
    logdensity_summand[l,] = dnbinom(l - 1, r, psi, log = T) + dcobin(x, theta, l, log = T)
  }
  logdensity = matrixStats::colLogSumExps(logdensity_summand)
  if(log){
    return(logdensity)
  }else{
    return(exp(logdensity))
  }
}


pmicobin <- function(q, theta, psi, r = 2, l_max = 70){
  if(length(psi) != 1){
    stop("psi must be scalar")
  }
  if(length(theta) != 1){
    stop("theta must be scalar")
  }
  logcdf_summand = matrix(-Inf, l_max, length(q))
  for(l in 1:l_max){
    logcdf_summand[l,] = dnbinom(l - 1, r, psi, log = T) + log(pcobin(q, theta, l))
  }
  logcdf = matrixStats::colLogSumExps(logcdf_summand)
  return(exp(logcdf))
}

# xgrid = seq(0.0001, 0.9999, length.out = 200)
# theta  = 5
# psi = 0.2
# draws = rcobin(10000, theta, 5)
# hist(draws)
# plot(ecdf(draws))
# lines(xgrid, pcobin(xgrid, theta, 5), col = 2); abline(v=0); abline(v=1)
# 
# draws = rmicobin(10000, theta, psi)
# hist(draws)
# plot(ecdf(draws))
# lines(xgrid, pmicobin(xgrid, theta, psi), col = 2); abline(v=0); abline(v=1)

emicobin <- function(eta){
  bftprime(eta)
}

vmicobin <- function(eta, psi){
  bftprimeprime(eta)*psi
}

