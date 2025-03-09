#' Vectorized sampling (by row by row) with probability stored in rows of matrix
#' source: https://stackoverflow.com/questions/20508658/sampling-repeatedly-with-different-probability
#'
#' @param Weight n by K matrix, each row corresponds to unnormalized probabilities of sampling weights
#'
#' @return integers with length n, each element in {1,...,K}
#' @export
#'
#' @examples
sample.rowwise <- function(Weight) {
  x <- runif(nrow(Weight))
  cumul.w <- Weight %*% upper.tri(diag(ncol(Weight)), diag = TRUE)/rowSums(Weight)
  i <- rowSums(x > cumul.w) + 1L
  i
}

# generalized logit transform, transform (a,b) range to real line
glogit <- function(x, xmin, xmax){
  x01 = (x - xmin)/(xmax - xmin)
  return(log(x01/(1-x01)))
}
# inverse of glogit
inv_glogit <- function(x, xmin, xmax){
  x01 = 1/(1+exp(-x))
  return(x01*(xmax - xmin) + xmin)
}


# ramdom sampling from KG(b,c) using infinite sum expression, truncated at kmax
rkg.gamma <- function(n, b= 1, z = 0, kmax = 200){
  k = 0:kmax # start from 0
  denom = (1+k)^2 + z^2/(4*pi^2) # length kmax + 1
  temp = matrix(rgamma(n*(kmax + 1), b, 1)/rep(denom, n), nrow = kmax + 1, ncol = n, byrow = F)
  colSums(temp)*(1/(2*pi^2))
}

