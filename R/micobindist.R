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
  logconst = pnbinom(l_max-1, r, psi, log = T)
  maxerror = (1-pnbinom(l_max-1, r, psi, log = FALSE))
  if(any(maxerror > 0.0001)) warning("psi is small, so that deviation from truncated and untruncated may be large")
  for(l in 1:l_max){
    logdensity_summand[l,] = dnbinom(l - 1, r, psi, log = T) - logconst + dcobin(x, theta, l, log = T)
  }
  logdensity = matrixStats::colLogSumExps(logdensity_summand)
  if(log){
    return(logdensity)
  }else{
    return(exp(logdensity))
  }
}
# 
# dmicobinold <- function(x, theta, psi, r = 2, log = FALSE, l_max = 70){
#   n = length(x)
#   if(length(psi) != 1){
#     stop("psi must be scalar")
#   }
#   if(length(theta) == 1) theta = rep(theta, length(x))
#   #if(length(psi) == 1) psi = rep(psi, length(x))
#   stopifnot("length of theta should be either 1 or length(x)" = (length(theta) == length(x)))
#   #stopifnot("length of psi should be either 1 or length(x)" = (length(psi) == length(x)))
#   
#   logdensity_summand = matrix(-Inf, l_max, length(x))
#   error = (1-pnbinom(l_max-1, r, psi, log = FALSE))
#   if(error > 0.0001) warning("m_max is too small and density may be off by 0.0001. Increase l_max.")
#   for(l in 1:l_max){
#     logdensity_summand[l,] = dnbinom(l - 1, r, psi, log = T) + dcobin(x, theta, l, log = T)
#   }
#   logdensity = matrixStats::colLogSumExps(logdensity_summand)
#   if(log){
#     return(logdensity)
#   }else{
#     return(exp(logdensity))
#   }
# }
# 
# dcobin(c(0.5,0.2), theta = c(2,2), lambda = c(2,2))
# xgrid = seq(0, 1, length = 1000)
# plot(xgrid, dmicobin(xgrid, theta = 2, 1/2, log = F))
# lines(xgrid, dmicobinold(xgrid, theta = 2, 1/2, log = F), col = 2)


#' Title
#'
#' @param q 
#' @param theta 
#' @param psi 
#' @param r 
#' @param l_max 
#'
#' @returns
#' @export
#'
#' @examples
pmicobin <- function(q, theta, psi, r = 2, l_max = 70){
  if(length(psi) != 1){
    stop("psi must be scalar")
  }
  if(length(theta) != 1){
    stop("theta must be scalar")
  }
  cdf_summand = matrix(0, l_max, length(q))
  error = (1-pnbinom(l_max-1, r, psi, log = FALSE))
  if(error > 0.0001) warning("l_max may be too small")
  for(l in 1:l_max){
    cdf_summand[l,] = dnbinom(l - 1, r, psi)*pcobin(q, theta, l)
  }
  cdf = colSums(cdf_summand)
  return(cdf)
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

