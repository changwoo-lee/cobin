# Vectorized sampling (by row by row) with probability stored in rows of matrix
# source: https://stackoverflow.com/questions/20508658/sampling-repeatedly-with-different-probability
#
# n by K matrix, each row corresponds to unnormalized probabilities of sampling weights
#
# return integers with length n, each element in {1,...,K}
sample.rowwise <- function(Weight) {
  W   <- as.matrix(Weight)
  rs  <- rowSums(W)
  if (any(!is.finite(rs) | rs <= 0)) stop("Row sums must be positive/finite.")
  P   <- W / rs
  CDF <- rowCumsums(P)
  CDF[, ncol(CDF)] <- 1            
  u   <- runif(nrow(CDF))          
  pmin.int(rowSums(u > CDF) + 1L, ncol(CDF))  # clamp
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
# 
# ekg = function(b, c){
#   if(c==0){
#     return( b/12 )
#   }else{
#     return( b*((c/2)*(1/tanh(c/2))-1)/(c^2) )
#   }
# }

ekg <- function(b, c) {
  nb <- length(b); nc <- length(c)
  
  # Broadcast scalars
  if (nb == 1L && nc > 1L) b <- rep(b, nc)
  if (nc == 1L && nb > 1L) c <- rep(c, nb)
  
  # Require matching lengths after broadcasting
  if (length(b) != length(c)) {
    stop("Lengths of b and c must match, or one of them must be length 1.")
  }
  
  out  <- numeric(length(c))
  zero <- (c == 0)
  
  # c == 0 branch: limit is b/12
  out[zero] <- b[zero] / 12
  
  # c != 0 branch
  cz <- c[!zero]
  out[!zero] <- b[!zero] * ((cz/2) * (1 / tanh(cz/2)) - 1) / (cz^2)
  
  out
}




