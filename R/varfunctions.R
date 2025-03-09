#' Cumulant (log partition) function of cobin
#'
#' B(x) = log((exp(x)-1)/x)
#'
#' @param x
#'
#' @returns
#' @export
#'
#' @examples
bft <- function(x){
  # divided domain into 5 different parts to maintain numerical stability
  zeroidx = which(abs(x) <= 1e-16) # practically zero
  pos1idx = which(x > 1e-16 & x < 1)
  pos2idx = which(x >= 1)
  neg1idx = which(x < -1e-16 & x > -1)
  neg2idx = which(x <= -1)
  out = numeric(length(x))
  if(length(zeroidx) > 0) out[zeroidx] = 0
  if(length(neg1idx) > 0) out[neg1idx] = log(-x[neg1idx]) - log(-expm1(x[neg1idx]))
  if(length(neg2idx) > 0) out[neg2idx] = log(x[neg2idx]/(exp(x[neg2idx])-1))
  if(length(pos1idx) > 0) out[pos1idx] = log(x[pos1idx]) - log(expm1(x[pos1idx]))
  if(length(pos2idx) > 0) out[pos2idx] = log(x[pos2idx]) - x[pos2idx]-log1p(-exp(-x[pos2idx]))
  return(-out) # return with opposite sign; above is log(x/(exp(x)-1))
}

#' First derivative of cobin cumulant (log partition) function
#'
#' B'(x) = 1/(1-exp(-x))-1/x.
#' When g is canonical link, this is same as g^(-1)
#'
#' @param x
#'
#' @returns
#' @export
#'
#' @examples
bftprime = function(x){
  #out = 1/(1-exp(-x))-1/x
  out = -1/expm1(-x)-1/x
  # if any element is NaN, replace NaN to 0.5
  out[is.nan(out)] = 0.5
  out
}

#' Second derivative of cobin cumulant (log partition) function
#'
#' B''(x) = 1/x^2 + 1/(2-2*cosh(x))
#' used Taylor series expansion for x near 0 for stability
#' When g is canonical link, this is same as (g^(-1))'
#'
#' @param x
#'
#' @returns
#' @export
#'
#' @examples
bftprimeprime = function(x){
  #out = 1/x^2 + 1/(2-2*cosh(x))
  out = 1/x^2 - 1/(expm1(x) + expm1(-x))
  nearzeroidx = which(abs(x) < 5e-3)
  out[nearzeroidx] = 1/12 - x[nearzeroidx]^2/240 + x[nearzeroidx]^4/6048
  out
}


#' Third derivative of cobin cumulant (log partition) function
#'
#' B'''(x) = 1/(4*tanh(x/2)*sinh(x/2)^2) - 2/x^3
#' used Taylor series expansion for x near 0 for stability
#'
#' @param x
#'
#' @returns
#' @export
#'
#' @examples
bftprimeprimeprime = function(x){
  out = 1/(4*tanh(x/2)*sinh(x/2)^2) - 2/x^3
  nearzeroidx = which(abs(x) < 5e-3)
  out[nearzeroidx] = -x[nearzeroidx]/120 + x[nearzeroidx]^3/1512 - x[nearzeroidx]^5/28800
  out
}

#' Inverse of first derivative of cobin cumulant (log partition) function
#'
#' Calculates (B')^(-1)(y) using numerical inversion (Newton-Raphson)
#' This is the cobit link function g, the canonical link function of cobin
#'
#' @param y
#' @param x0 initial value
#' @param tol tolerance, stopping criterion for Newton-Raphson
#' @param maxiter max iteration of Newton-Raphson
#'
#' @returns
#' @export
#'
#' @examples
bftprimeinv = function(y, x0 = 0, tol = 1e-8, maxiter = 100){
  # y is vector
  # check if all elements of y is between 0 and 1
  if(any(y <= 0) || any(y >= 1)){
    stop("y must be between 0 and 1")
  }
  x = rep(x0, length(y))
  xnew = rep(NA, length(y))
  update =  rep(T, length(y))
  for(i in 1:maxiter){
    idx = which(update)
    xnew[idx] = x[idx] - (bftprime(x[idx])-y[idx])/bftprimeprime(x[idx])
    update = abs(x - xnew) > tol
    x = xnew
    if(sum(update)==0){ # all FALSE
      break
    }
  }
  return(x)
}


#' Variance function of cobin
#'
#' B''(B'^(-1)(mu))
#'
#' @param mu
#'
#' @returns
#' @export
#'
#' @examples
Vft <- function(mu){
  bftprimeprime(bftprimeinv(mu))
}


# gprime <- function(y){
# #   1/ginvprime(g(y))
#   1/bftprimeprime(bftprimeinv(y))
# }





#
# # log(x/(exp(x)-1)),
#
# log_x_over_expxm1 <- function(x){
#   # divided domain into 5 different parts to maintain numerical stability
#   zeroidx = which(abs(x) <= 1e-16) # practically zero
#   pos1idx = which(x > 1e-16 & x < 1)
#   pos2idx = which(x >= 1)
#   neg1idx = which(x < -1e-16 & x > -1)
#   neg2idx = which(x <= -1)
#   out = numeric(length(x))
#   if(length(zeroidx) > 0) out[zeroidx] = 0
#   if(length(neg1idx) > 0) out[neg1idx] = log(-x[neg1idx]) - log(-expm1(x[neg1idx]))
#   if(length(neg2idx) > 0) out[neg2idx] = log(x[neg2idx]/(exp(x[neg2idx])-1))
#   if(length(pos1idx) > 0) out[pos1idx] = log(x[pos1idx]) - log(expm1(x[pos1idx]))
#   if(length(pos2idx) > 0) out[pos2idx] = log(x[pos2idx]) - x[pos2idx]-log1p(-exp(-x[pos2idx]))
#   out
# }
# #xgrid = c(-3000, -2, -1e-13, 0, 2e-11, 4, 6000)
# #log(xgrid/(exp(xgrid)-1))
# #log_x_over_expxm1(xgrid)


