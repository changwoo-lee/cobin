
#' Random variate generation for continuous binomial distribution
#'
#' y ~ cobin(theta, 1/lambda)
#' which is distributionally equivalent to
#' x_1 + ... + x_lambda, x_l ~ contiBernoulli(theta)
#'
#' Continuous Bernoulli distribution with parameter theta has a density function
#' f(y; theta) = theta/(e^theta-1) * e^(theta * y)
#' or equivalently, with phi = e^theta/(1+e^theta),
#' f(y; phi) \propto phi^y * (1-phi)^(1-y)
#'
#' @param n integer, number of samples
#' @param theta scalar or length n vector, theta
#' @param lambda scalar or length n vector, lambda, length should be same as theta
#'
#' @returns
#' @export
#'
#' @examples
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










dcobin <- function(x, theta, lambda, log = FALSE){
  n = length(x)
  if(length(theta) == 1){
    theta = rep(theta, n)
  }
  # length of lambda must be 1
  if(length(theta) != n | length(lambda) != 1){
    stop("length of theta should be either 1 or length of x, and leegnth of lambda must be 1")
  }
  logdensity = log(lambda) + dIH(lambda*x, lambda, log = T) + lambda*theta*x - lambda*bft(theta)
  if(log){
    return(logdensity)
  }else{
    return(exp(logdensity))
  }
}






ecobin <- function(theta){
  bftprime(theta)
}

vcobin <- function(theta, lambda){
  bftprimeprime(theta)/lambda
}
