#' Random variate generation for cobin (continuous binomial) distribution
#'
#' Continuous binomial distribution with natural parameter \eqn{\theta} and dispersion parameter \eqn{1/\lambda}, in short \eqn{Y \sim cobin(\theta, \lambda^{-1})}, has density 
#' \deqn{
#'  p(y; \theta, \lambda^{-1}) = h(y;\lambda) \exp(\lambda \theta y - \lambda B(\theta)), \quad 0 \le y \le 1
#' }
#' where \eqn{B(\theta) = \log\{(e^\theta - 1)/\theta\}} and \eqn{h(y;\lambda) = \frac{\lambda}{(\lambda-1)!}\sum_{k=0}^{\lambda} (-1)^k {\lambda \choose k} \max(0,\lambda y-k)^{\lambda-1}}. 
#' When \eqn{\lambda = 1}, it becomes continuous Bernoulli distribution. 
#' 
#' The random variate generation is based on the fact that \eqn{cobin(\theta, \lambda^{-1})} is equal in distribution to the sum of \eqn{\lambda} \eqn{cobin(\theta, 1)} random variables, scaled by \eqn{\lambda^{-1}}.
#' Random variate generation for continuous Bernoulli is done by inverse cdf transform method. 
#'
#' @param n integer, number of samples
#' @param theta scalar or length n vector, natural parameter.
#' @param lambda scalar or length n vector, inverse of dispersion parameter. Must be integer, length should be same as theta
#'
#' @returns random samples from \eqn{cobin(\theta,\lambda^{-1})}.
#' @export
#'
#' @examples
#' hist(rcobin(1000, 2, 3), freq = FALSE)
#' xgrid = seq(0, 1, length = 500)
#' lines(xgrid, dcobin(xgrid, 2, 3))
#' 
rcobin <- function(n, theta, lambda){
  # check lambda is integer
  if(any(lambda %% 1 != 0)){
    stop("lambda must be integer")
  }
  if(length(theta)==1 & length(lambda)==1){
    return(colMeans(matrix(rcb(n*lambda, theta), lambda, n)))
  }else if(length(theta) == n & length(lambda) == n){
    rep_idx = rep.int(seq_len(n), lambda) # if m = c(3,1,2), this creates c(1,1,1,2,3,3)
    draws = rcb(sum(lambda), theta[rep_idx])
    #return(Rfast::group(draws, rep_idx, method = "mean"))
    #return(as.numeric(rowsum(as.matrix(draws), group = rep_idx)/tabulate(rep_idx)))
    return(as.numeric(rowsum(as.matrix(draws), group = rep_idx)/lambda))
  }else{
    stop("length of theta and lambda must be both 1 or n")
  }
}










#' Density function of cobin (continuous binomial) distribution
#' 
#' Continuous binomial distribution with natural parameter \eqn{\theta} and dispersion parameter \eqn{1/\lambda}, in short \eqn{Y \sim cobin(\theta, \lambda^{-1})}, has density 
#' \deqn{
#'  p(y; \theta, \lambda^{-1}) = h(y;\lambda) \exp(\lambda \theta y - \lambda B(\theta)), \quad 0 \le y \le 1
#' }
#' where \eqn{B(\theta) = \log\{(e^\theta - 1)/\theta\}} and \eqn{h(y;\lambda) = \frac{\lambda}{(\lambda-1)!}\sum_{k=0}^{\lambda} (-1)^k {\lambda \choose k} \max(0,\lambda y-k)^{\lambda-1}}. 
#' When \eqn{\lambda = 1}, it becomes continuous Bernoulli distribution. 
#' 
#' For the evaluation of \eqn{h(y;\lambda)}, see `?cobin::dIH`. 
#' 
#' @param x num (length n), between 0 and 1, evaluation point
#' @param theta scalar or length n vector, num (length 1 or n), natural parameter
#' @param lambda scalar or length n vector, integer, inverse of dispersion parameter
#' @param log logical (Default FALSE), if TRUE, return log density
#'
#' @returns density of \eqn{cobin(\theta,\lambda^{-1})}
#' @export
#'
#'
#' @examples
#' 
#' xgrid = seq(0, 1, length = 500)
#' plot(xgrid, dcobin(xgrid, 0, 1), type="l", ylim = c(0,3)) # uniform 
#' lines(xgrid, dcobin(xgrid, 0, 3))
#' plot(xgrid, dcobin(xgrid, 2, 3), type="l")
#' lines(xgrid, dcobin(xgrid, -2, 3))
#' 
dcobin <- function(x, theta, lambda, log = FALSE){
  n = length(x)
  if(length(theta) == 1) theta = rep(theta, n)
  if(length(lambda) == 1) lambda = rep(lambda, n)
  if(any(lambda > 80)){
    warning("density calculation with lambda > 80 is unstable")
  }
  stopifnot("length of theta should be either 1 or length(x)" = (length(theta) == length(x)))
  stopifnot("length of lambda should be either 1 or length(x)" = (length(lambda) == length(x)))
  # # length of lambda must be 1
  # if(length(theta) != n | length(lambda) != 1){
  #   stop("length of theta should be either 1 or length of x, and leegnth of lambda must be 1")
  # }
  logdensity = log(lambda) + dIH(lambda*x, lambda, log = TRUE) + lambda*theta*x - lambda*bft(theta)
  if(log){
    return(logdensity)
  }else{
    return(exp(logdensity))
  }
}



#' Cumulative distribution function of cobin (continuous binomial) distribution
#' 
#' Continuous binomial distribution with natural parameter \eqn{\theta} and dispersion parameter \eqn{1/\lambda}, in short \eqn{Y \sim cobin(\theta, \lambda^{-1})}, has density 
#' \deqn{
#'  p(y; \theta, \lambda^{-1}) = h(y;\lambda) \exp(\lambda \theta y - \lambda B(\theta)), \quad 0 \le y \le 1
#' }
#' where \eqn{B(\theta) = \log\{(e^\theta - 1)/\theta\}} and \eqn{h(y;\lambda) = \frac{\lambda}{(\lambda-1)!}\sum_{k=0}^{\lambda} (-1)^k {\lambda \choose k} \max(0,\lambda y-k)^{\lambda-1}}. 
#' When \eqn{\lambda = 1}, it becomes continuous Bernoulli distribution. 
#' 
#' @param q num (length n), between 0 and 1, evaluation point
#' @param theta scalar, natural parameter
#' @param lambda integer, inverse of dispersion parameter
#'
#' @returns c.d.f. of \eqn{cobin(\theta,\lambda^{-1})}
#' @export
#'
#' @examples
#' 
#' xgrid = seq(0, 1, length = 500)
#' out = pcobin(xgrid, 1, 2)
#' plot(ecdf(rcobin(10000, 1, 2)))
#' lines(xgrid, out, col = 2)
#' 
pcobin <- function(q, theta, lambda){
  if(length(lambda)!=1) stop("lambda must be scalar")
  # check lambda is integer
  if(lambda %% 1 != 0){
    stop("lambda must be integer")
  }
  if(lambda == 1) return(pcb(q, theta))
  if(length(theta)!=1) stop("theta must be scalar")
  # check x is between 0 and 1
  if(any(q < 0 | q > 1)){
    stop("q must be between 0 and 1")
  }
  # numerical integration for lambda > 40
  if(lambda > 40 & lambda <= 60){
    return(pcobin_numerical(q, theta, lambda))
  }else if(lambda > 60){
    #warning("normal cdf approximation used for lambda > 60 due to numerical stability")
    return(pnorm(q, mean = bftprime(theta), sd = sqrt(vcobin(theta,lambda))))
  }
  if(theta < 0){
    out = pcobin_negtheta(q, theta, lambda)
  }else if(theta == 0){
    # to imporve numerical stability
    idx1 = which(q <= 0.5)
    idx2 = which(q > 0.5)
    out = rep(0, length(q))
    if(length(idx1) > 0) out[idx1] = pIH01(q[idx1], lambda)
    if(length(idx2) > 0) out[idx2] = 1 - pIH01(1-q[idx2], lambda)
  }else{
    out = 1-pcobin_negtheta(1-q, -theta, lambda)
  }
  return(pmin(pmax(0, out), 1))
}

pIH01 <- function(q, lambda, log = F){
  n = length(q)
  logsummand_mat = matrix(0, lambda+1, n)
  common = log(lambda) - lfactorial(lambda-1) + lchoose(lambda, 0:lambda) 
  for(i in 1:n){
    temp = lambda*q[i] - (0:lambda)
    zeroidx = which(temp <= 0)
    nonzeroidx = which(temp > 0)
    logsummand_mat[zeroidx,i] = -Inf
    logsummand_mat[nonzeroidx,i] = common[nonzeroidx] + lambda*log(temp[nonzeroidx]) - 2*log(lambda)
  }
  signs = (-1)^(0:lambda)
  logsums_positive = matrixStats::colLogSumExps(logsummand_mat[signs == 1,, drop = F])
  logsums_negative = matrixStats::colLogSumExps(logsummand_mat[signs == -1,, drop = F])
  if(any(logsums_positive < logsums_negative)){
    warning("numerical error, return 0 density value")
    logsums_positive = logsums_negative + pmax(logsums_positive-logsums_negative, 0)
    
    logcdf = logsums_positive + log1p(-exp(logsums_negative-logsums_positive)) # log-minus-exp
  }else{
    logcdf = logsums_positive + log1p(-exp(logsums_negative-logsums_positive)) # log-minus-exp
  }
  idxnan = which(is.nan(logcdf))
  logcdf[idxnan] = -Inf
  if(log){
    return(logcdf)
  }else{
    return(exp(logcdf))
  }
}

pcobin_negtheta2 <- function(q, theta, lambda, log = FALSE){
  n = length(q)
  summand_mat = matrix(0, lambda+1, n)
  common = log(lambda) + lchoose(lambda, 0:lambda) - lambda*bft(theta)
  for(i in 1:n){
    temp = lambda*q[i] - (0:lambda)
    zeroidx = which(temp <= 0)
    nonzeroidx = which(temp > 0)
    summand_mat[zeroidx,i] = 0
    # https://dl.acm.org/doi/10.1145/2972951#supplementary-materials
    # lower incomplete gamma function in R is only supported by positive argument; (unless using gsl::hyperg_1F1 which is numerically instable)
    # currently, lambda > 40 cause numerical issues. we can apply same trick of pIH (symmetry) to avoid this issue if lower incomplete gamma can be evaluated in negative domain
    #term1 = pgamma( - theta * (lambda * q[i] - (0:lambda)[nonzeroidx]), lambda)#*gamma(lambda) # only valid when theta < 0 
    logterm1 = pgamma( - theta * (lambda * q[i] - (0:lambda)[nonzeroidx]), lambda, log.p = TRUE)#*gamma(lambda) # only valid when theta < 0 
    
    #term2 = pgamma( - theta * (lambda * (0:lambda)[nonzeroidx]/lambda - (0:lambda)[nonzeroidx]), lambda)#*gamma(lambda) # only valid when theta < 0 
    #summand_mat[nonzeroidx,i] = exp(common[nonzeroidx])*exp(theta*(0:lambda)[nonzeroidx])/(lambda)*term1/abs(theta)^lambda
    summand_mat[nonzeroidx,i] = exp(common[nonzeroidx] + theta*(0:lambda)[nonzeroidx] - log(lambda) + logterm1 - lambda*log(abs(theta)))
  }
  signs = (-1)^(0:lambda)
  #browser()
  return(colSums(summand_mat*signs))
}

pcobin_negtheta <- function(q, theta, lambda, log = F){
  n = length(q)
  logsummand_mat = matrix(0, lambda+1, n)
  common = log(lambda) - lfactorial(lambda-1) + lchoose(lambda, 0:lambda) - lambda*bft(theta)
  for(i in 1:n){
    temp = lambda*q[i] - (0:lambda)
    zeroidx = which(temp <= 0)
    nonzeroidx = which(temp > 0)
    logsummand_mat[zeroidx,i] = -Inf
    # https://dl.acm.org/doi/10.1145/2972951#supplementary-materials
    # lower incomplete gamma function in R is only supported by positive argument; (unless using gsl::hyperg_1F1 which is numerically instable)
    # currently, lambda > 40 cause numerical issues. we can apply same trick of pIH (symmetry) to avoid this issue if lower incomplete gamma can be evaluated in negative domain
    
    #term1 = pgamma( - theta * (lambda * q[i] - (0:lambda)[nonzeroidx]), lambda)*gamma(lambda) # only valid when theta < 0 
    #term2 = pgamma( - theta * (lambda * (0:lambda)[nonzeroidx]/lambda - (0:lambda)[nonzeroidx]), lambda)*gamma(lambda) # only valid when theta < 0 
    #logsummand_mat[nonzeroidx,i] = common[nonzeroidx] + theta*(0:lambda)[nonzeroidx] - log(lambda) + log((term1 - term2)/abs(theta)^lambda)
    logterm1 = pgamma( - theta * (lambda * q[i] - (0:lambda)[nonzeroidx]), lambda, log.p = TRUE)#*gamma(lambda) # only valid when theta < 0 
    #term2 = pgamma( - theta * (lambda * (0:lambda)[nonzeroidx]/lambda - (0:lambda)[nonzeroidx]), lambda)#*gamma(lambda) # only valid when theta < 0 
    logsummand_mat[nonzeroidx,i] = common[nonzeroidx] + theta*(0:lambda)[nonzeroidx] - log(lambda) + lgamma(lambda) + logterm1 - lambda*log(abs(theta))
  }
  signs = (-1)^(0:lambda)
  logsums_positive = matrixStats::colLogSumExps(logsummand_mat[signs == 1,, drop = F])
  logsums_negative = matrixStats::colLogSumExps(logsummand_mat[signs == -1,, drop = F])
  if(any(logsums_positive < logsums_negative)){
    warning("numerical error, return 0 density value")
    logsums_positive = logsums_negative + pmax(logsums_positive-logsums_negative, 0)
    
    logcdf = logsums_positive + log1p(-exp(logsums_negative-logsums_positive)) # log-minus-exp
  }else{
    logcdf = logsums_positive + log1p(-exp(logsums_negative-logsums_positive)) # log-minus-exp
  }
  idxnan = which(is.nan(logcdf))
  logcdf[idxnan] = -Inf
  if(log){
    return(logcdf)
  }else{
    return(exp(logcdf))
  }
}

pcobin_numerical <- function(q, theta, lambda) {
  sapply(q, function(q_val) {
    integrate(dcobin, 0, q_val, theta, lambda, subdivisions = 1000L)$value
  })
}

#check
# q = seq(0.0001, 0.9999, length = 5)
# theta = 1
# lambda = 50
# pcobin(q, theta, lambda)
# pcobin_hypergeometric(q, theta, lambda)
# pcobin_numerical(q, theta, lambda)
# 
ecobin <- function(theta){
  bftprime(theta)
}

vcobin <- function(theta, lambda){
  bftprimeprime(theta)/lambda
}

