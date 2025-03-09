#' Density of Irwin-hall distribution
#'
#' Irwin-Hall distribution is
#'
# avoid overflow by using log scale
# and using symmetry at x = m/2
# WARNING: m > 70~80 might have numerical issues\
#'
#' @param x
#' @param m
#' @param log
#'
#' @returns
#' @export
#'
#' @examples
dIH <- function(x, m, log = F){
  n = length(x)
  # check x is between 0 and 1
  if(any(x < 0 | x > m)){
    stop("x must be between 0 and m")
  }
  if(m==1){
    if(log) return(rep(0,n)) else return(rep(1,n))
  }
  idx = which(x > m/2)
  x[idx] = m - x[idx]
  logsummand_mat = matrix(0, m+1, n)
  temp = -lfactorial(m-1) + lchoose(m, 0:m)
  for(i in 1:n){
    logsummand_mat[,i] = temp + (m-1)*log(pmax(x[i]-0:m, 0))
  }
  signs = (-1)^(0:m)
  logsums_positive = matrixStats::colLogSumExps(logsummand_mat[signs == 1,, drop = F])
  logsums_negative = matrixStats::colLogSumExps(logsummand_mat[signs == -1,, drop = F])
  if(any(logsums_positive < logsums_negative)){
    warning("numerical error, return 0 density value")
    logsums_positive = logsums_negative + pmax(logsums_positive-logsums_negative, 0)

    logdensity = logsums_positive + log1p(-exp(logsums_negative-logsums_positive)) # log-minus-exp
  }else{
    logdensity = logsums_positive + log1p(-exp(logsums_negative-logsums_positive)) # log-minus-exp
  }
  idxnan = which(is.nan(logdensity))
  logdensity[idxnan] = -Inf
  if(log){
    return(logdensity)
  }else{
    return(exp(logdensity))
  }
}
#
# dIH(x = 0.3, m = 96, log = F)
# # m = 60
# # xgrid= seq(0, m, length = 10000)
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

