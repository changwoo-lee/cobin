#' Cumulant (log partition) function of cobin
#'
#' \eqn{B(x) = \log\{(\exp(x)-1)/x)\}}
#'
#' @param x input vector
#'
#' @returns \eqn{B(x) = \log\{(\exp(x)-1)/x)\}}
#' @export
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
#' \eqn{B'(x) = 1/(1-\exp(-x))-1/x}.
#' When \eqn{g} is canonical link of cobin, this is same as \eqn{g^{-1}}
#'
#' @param x input vector
#'
#' @returns \eqn{B'(x) = 1/(1-\exp(-x))-1/x}.
#' @export
bftprime = function(x){
  #out = 1/(1-exp(-x))-1/x
  out = -1/expm1(-x)-1/x
  # if any element is NaN, replace NaN to 0.5
  out[is.nan(out)] = 0.5
  out
}

#' @rdname bftprime
#' @export
cobitlinkinv <- bftprime

#' Second derivative of cobin cumulant (log partition) function
#'
#' \eqn{B''(x) = 1/x^2 + 1/(2-2*\cosh(x))}
#' used Taylor series expansion for x near 0 for stability
#'
#' @param x input vector
#'
#' @returns \eqn{B''(x) = 1/x^2 + 1/(2-2*\cosh(x))}
#' @export
bftprimeprime = function(x){
  #out = 1/x^2 + 1/(2-2*cosh(x))
  out = 1/x^2 - 1/(expm1(x) + expm1(-x))
  nearzeroidx = which(abs(x) < 5e-3)
  out[nearzeroidx] = 1/12 - x[nearzeroidx]^2/240 + x[nearzeroidx]^4/6048
  out
}


#' Third derivative of cobin cumulant (log partition) function
#'
#' \eqn{B'''(x) = 1/(4*\tanh(x/2)*\sinh(x/2)^2) - 2/x^3}
#' used Taylor series expansion for x near 0 for stability
#'
#' @param x input vector
#'
#' @returns \eqn{B'''(x) = 1/(4*\tanh(x/2)*\sinh(x/2)^2) - 2/x^3}
#' @export
bftprimeprimeprime = function(x){
  out = 1/(4*tanh(x/2)*sinh(x/2)^2) - 2/x^3
  nearzeroidx = which(abs(x) < 5e-3)
  out[nearzeroidx] = -x[nearzeroidx]/120 + x[nearzeroidx]^3/1512 - x[nearzeroidx]^5/28800
  out
}

#' Inverse of first derivative of cobin cumulant (log partition) function
#'
#' Calculates \eqn{(B')^{-1}(y)} using numerical inversion (Newton-Raphson),
#' where \eqn{B'(x) = 1/(1-\exp(-x))-1/x}. 
#' This is the cobit link function g, the canonical link function of cobin. 
#'
#' @param y input vector
#' @param x0 Defult 0, initial value
#' @param tol tolerance, stopping criterion for Newton-Raphson
#' @param maxiter max iteration of Newton-Raphson
#'
#' @returns \eqn{(B')^{-1}(y)}
#' @export
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

#' @rdname bftprimeinv
#' @export
cobitlink <- bftprimeinv

#' Variance function of cobin
#'
#' \eqn{B''(B'^{-1}(\mu))}
#'
#' @param mu input vector
#'
#' @returns \eqn{B''(B'^{-1}(\mu))}
#' @export
Vft <- function(mu){
  bftprimeprime(bftprimeinv(mu))
}



