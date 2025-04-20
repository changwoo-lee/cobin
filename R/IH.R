#' Density of Irwin-Hall distribution
#'
#' Irwin-Hall distribution is defined as a sum of m uniform (0,1) distribution. 
#' Its density is given as 
#' \deqn{
#'   f(x;m) = \frac{1}{(m-1)!}\sum_{k=0}^{m} (-1)^k {m \choose k} \max(0,x-k)^{m-1}, 0 < x < m
#' }
#' The density of Bates distribution, defined as an average of m uniform (0,1) distribution, can be obtained from change-of-variable (y = x/m),
#' \deqn{
#'   h(y;m) = \frac{m}{(m-1)!}\sum_{k=0}^{m} (-1)^k {m \choose k} \max(0,my-k)^{m-1}, 0 < y < 1
#' }
#' 
#' Due to alternating series representation, m > 80 may yield numerical issues
#'
#' @param x vector of quantities, between 0 and m
#' @param m integer, parameter
#' @param log logical, return log density if TRUE
#'
#' @returns (log) density evaluated at x
#' @importFrom matrixStats colLogSumExps
#' @export
#'
#' @examples
#' m = 4
#' xgrid= seq(0, m, length = 500)
#' plot(xgrid, dIH(xgrid, m, log = FALSE))
#' 
dIH <- function(x, m, log = F){
  n = length(x)
  if(length(m)==1) m = rep(m, n)
  stopifnot("length of m should be either 1 or length(x)" = (length(m) == length(x)))
  stopifnot(all(x >= 0 & x <= m))
  # check x is between 0 and 1
  logout = rep(0,n)
  notoneidx = which(m>1)
  if(length(notoneidx)==0){
    if(log) return(logout) else return(exp(logout))
  }
  # dealing with m > 1
  x = x[notoneidx]
  m = m[notoneidx]
  n = length(x)
  mmax = max(m)
  x = pmin(x, m-x)
  xmat = matrix(x, mmax + 1, n, byrow = TRUE)
  mmat = matrix(m, mmax + 1, n, byrow = TRUE)
  mseqmat = matrix(0:mmax, mmax + 1, n)
  #logsummand_mat[(1:(m[i]+1)),i] = -lfactorial(m[i]-1) + lchoose(m[i], 0:m[i]) + (m[i]-1)*log(pmax(x[i]-0:m[i], 0))
  #logsummand_mat = -lfactorial(mmat - 1) + lchoose(mmat, mseqmat) + (mmat - 1)*log(pmax(xmat - mseqmat, 0))
  logsummand_mat = -matrix(lfactorial(m-1), mmax + 1, n, byrow = TRUE) + lchoose(mmat, mseqmat) + (mmat - 1)*log(pmax(xmat - mseqmat, 0))
  #browser()
  signs = (-1)^(0:mmax)
  logsums_positive = matrixStats::colLogSumExps(logsummand_mat[signs == 1,, drop = F])
  logsums_negative = matrixStats::colLogSumExps(logsummand_mat[signs == -1,, drop = F])
  if(any(logsums_positive < logsums_negative)){
    warning("numerical error, return 0 density value")
    logsums_positive = logsums_negative + pmax(logsums_positive-logsums_negative, 0)

    logdensity = logsums_positive + log1p(-exp(logsums_negative-logsums_positive)) # log-minus-exp
  }else{
    logdensity = logsums_positive + log1p(-exp(logsums_negative-logsums_positive)) # log-minus-exp
  }
  logout[notoneidx] = logdensity
  idxnan = which(is.nan(logout))
  logout[idxnan] = -Inf
  if(log) return(logout) else return(exp(logout))
}


