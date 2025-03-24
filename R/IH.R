#' Density of Irwin-Hall distribution
#'
#' Irwin-Hall distribution with parameter m is defined as a sum of m uniform (0,1) distribution. 
#' 
#' avoid overflow by using log scale
#' and using symmetry at x = m/2
#' WARNING: m > 70~80 might have numerical issues\
#'
#' @param x num(length n), between 0 and m, evaluation point
#' @param m integer, parameter
#' @param log logical, return log density if TRUE
#'
#' @returns density value at x
#' @importFrom matrixStats colLogSumExps
#' @export
#'
#' @examples
#' 
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

# 
# dIHold2 <- function(x, m, log = F){
#   n = length(x)
#   # check x is between 0 and 1
#   if(any(x < 0 | x > m)){
#     stop("x must be between 0 and m")
#   }
#   if(length(m)==1) m = rep(m, n)
#   logout = rep(0,n)
#   notoneidx = which(m>1)
#   if(length(notoneidx)==0){
#     if(log) return(logout) else return(exp(logout))
#   }
#   # dealing with m > 1
#   x = x[notoneidx]
#   m = m[notoneidx]
#   n = length(x)
#   mmax = max(m)
#   x = pmin(x, m-x)
#   logsummand_mat = matrix(-Inf, mmax+1, n)
#   for(i in 1:n){
#     logsummand_mat[(1:(m[i]+1)),i] = -lfactorial(m[i]-1) + lchoose(m[i], 0:m[i]) + (m[i]-1)*log(pmax(x[i]-0:m[i], 0))
#   }
#   signs = (-1)^(0:mmax)
#   logsums_positive = matrixStats::colLogSumExps(logsummand_mat[signs == 1,, drop = F])
#   logsums_negative = matrixStats::colLogSumExps(logsummand_mat[signs == -1,, drop = F])
#   if(any(logsums_positive < logsums_negative)){
#     warning("numerical error, return 0 density value")
#     logsums_positive = logsums_negative + pmax(logsums_positive-logsums_negative, 0)
#     
#     logdensity = logsums_positive + log1p(-exp(logsums_negative-logsums_positive)) # log-minus-exp
#   }else{
#     logdensity = logsums_positive + log1p(-exp(logsums_negative-logsums_positive)) # log-minus-exp
#   }
#   logout[notoneidx] = logdensity
#   idxnan = which(is.nan(logout))
#   logout[idxnan] = -Inf
#   if(log) return(logout) else return(exp(logout))
# }
# 
# # 
# dIHold <- function(x, m, log = F){
#   n = length(x)
#   # check x is between 0 and 1
#   if(any(x < 0 | x > m)){
#     stop("x must be between 0 and m")
#   }
#   if(m==1){
#     if(log) return(rep(0,n)) else return(rep(1,n))
#   }
#   idx = which(x > m/2)
#   x[idx] = m - x[idx]
#   logsummand_mat = matrix(0, m+1, n)
#   temp = -lfactorial(m-1) + lchoose(m, 0:m)
#   for(i in 1:n){
#     logsummand_mat[,i] = temp + (m-1)*log(pmax(x[i]-0:m, 0))
#   }
#   signs = (-1)^(0:m)
#   logsums_positive = matrixStats::colLogSumExps(logsummand_mat[signs == 1,, drop = F])
#   logsums_negative = matrixStats::colLogSumExps(logsummand_mat[signs == -1,, drop = F])
#   if(any(logsums_positive < logsums_negative)){
#     warning("numerical error, return 0 density value")
#     logsums_positive = logsums_negative + pmax(logsums_positive-logsums_negative, 0)
# 
#     logdensity = logsums_positive + log1p(-exp(logsums_negative-logsums_positive)) # log-minus-exp
#   }else{
#     logdensity = logsums_positive + log1p(-exp(logsums_negative-logsums_positive)) # log-minus-exp
#   }
#   idxnan = which(is.nan(logdensity))
#   logdensity[idxnan] = -Inf
#   if(log){
#     return(logdensity)
#   }else{
#     return(exp(logdensity))
#   }
# }

# 
# n = 5000
# x = runif(n, 0, 5)
# m = 50
# microbenchmark::microbenchmark(
# dIH(x, m, log = F),
# dIH2(x, m, log = F),
# dIHold(x, m, log = F))
# 
# mvec = sample(c(8, 9, 10), n, replace = T)
# 
# all.equal(dIH(x, mvec, log = F), dIH2(x, mvec, log = F))
# microbenchmark::microbenchmark(
# dIH(x, mvec, log = F),
# dIH2(x, mvec, log = F)
# )
# # # m = 60
# xgrid= seq(0, m, length = 10000)
# # plot(xgrid/m, dIH(xgrid, m, log = T))
# # m = 70
# # xgrid= seq(0, m, length = 10000)
# # lines(xgrid/m, dIH(xgrid, m, log = T), col = "red")
# # m = 80
# # xgrid= seq(0, m, length = 10000)
# # lines(xgrid/m, dIH(xgrid, m, log = T), col = "blue")
# m = 200
# xgrid= seq(0, m, length = 10000)
# plot(xgrid/m, dIH(xgrid/3, m, log = T), col = "green")

