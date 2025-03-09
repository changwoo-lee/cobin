
#' Quantile function of continuous Bernoulli (continuous binomial with lambda = 1)
#'
#' Continuous Bernoulli distribution with parameter theta has a density function
#' f(y; theta) = theta/(e^theta-1) * e^(theta * y)
#' or equivalently, with phi = e^theta/(1+e^theta),
#' f(y; phi) \propto phi^y * (1-phi)^(1-y)
#'
#' Its quantile function is
#'
#' log((exp(theta)-1)*p+1)/theta
#'
#' below codes evalualtes in a numerically stable way
#'
#' @param p length n vector, probabilities
#' @param theta scalar or length n vector, theta
#'
#' @returns length n vector of quantiles
#' @export
#'
#' @examples
#'
#' pgrid = seq(0.001, 0.999, length = 1000)
#' plot(pgrid, qcb(pgrid, 0)) # quantile function of uniform
#' plot(pgrid, qcb(pgrid, 3)) #
#'
qcb <- function(p, theta){
  if(length(p) != length(theta) & length(theta) != 1 ){
    stop("length of p and theta must be the same, or theta has length 1")
  }
  theta = array(theta,length(p))
  zeroidx = which(abs(theta) < 1e-12)
  posidx = which(theta > 1e-12)
  negidx = which(theta < -1e-12)
  out = numeric(length(p))
  out[zeroidx] = p[zeroidx]
  out[posidx] = 1 + log(p[posidx]+(1-p[posidx])*exp(-theta[posidx]))/theta[posidx]
  out[negidx] = log((exp(theta[negidx])-1)*p[negidx]+1)/theta[negidx]
  return(out)
}

# random variate generation of continuous Bernoulli (continuous binomial with lambda = 1)
rcb <- function(n, theta){
  return(qcb(p = runif(n), theta))
}
