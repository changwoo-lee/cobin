build_Q_exponential <- function(distmat, sigma.sq, phi,
                                 Nlist, ord){
  R = sigma.sq * exp(-distmat * phi) # exponential kernel
  n = ncol(R)
  R = R[ord,ord] # important! order
  # sparse matrix
  A = matrix(0, n, n)
  D = numeric(n) #D = Matrix(0, n, n);
  D[1] = 1
  # pseudocode 2 of Finley et al.
  #Nlist[[i]] be the set of indices j  <i such that A[i,j] not qe 0
  for(i in 1:(n-1)) {
    temp =  solve(R[Nlist[[i+1]], Nlist[[i+1]]],R[Nlist[[i+1]], i+1])
    A[i+1, Nlist[[i+1]]]= temp
    D[i+1]= R[i+1, i+1] - sum(R[i+1, Nlist[[i+1]]]*temp)
  }
  A = as(A, "dgCMatrix")

  Dmat = Matrix::Diagonal(n, x = D)
  Q_reordered = Matrix::t((Matrix::Diagonal(n) - A))%*%Matrix::solve(Dmat)%*%(Matrix::Diagonal(n) - A)
  #Linv =  Matrix::Diagonal(n, x = 1/sqrt(D))%*%(Matrix::Diagonal(n) - A)
  Q = Q_reordered[order(ord),order(ord)]
  Q = as(Q, "dgCMatrix")
  return(Q)
}
# 
# build_Q_exponential <- function(distmat, phi,
#                                 Nlist, ord){
#   R = exp(-distmat * phi) # exponential kernel
#   n = ncol(R)
#   R = R[ord,ord] # important! order
#   # sparse matrix
#   A = matrix(0, n, n)
#   D = numeric(n) #D = Matrix(0, n, n);
#   D[1] = 1
#   # pseudocode 2 of Finley et al.
#   #Nlist[[i]] be the set of indices j  <i such that A[i,j] not qe 0
#   for(i in 1:(n-1)) {
#     temp =  solve(R[Nlist[[i+1]], Nlist[[i+1]]],R[Nlist[[i+1]], i+1])
#     A[i+1, Nlist[[i+1]]]= temp
#     D[i+1]= R[i+1, i+1] - sum(R[i+1, Nlist[[i+1]]]*temp)
#   }
#   A = as(A, "dgCMatrix")
#   Dmat = Matrix::Diagonal(n, x = D)
#   Q_reordered = Matrix::t((Matrix::Diagonal(n) - A))%*%Matrix::solve(Dmat)%*%(Matrix::Diagonal(n) - A)
#   #Linv =  Matrix::Diagonal(n, x = 1/sqrt(D))%*%(Matrix::Diagonal(n) - A)
#   Q = Q_reordered[order(ord),order(ord)]
#   Q = as(Q, "dgCMatrix")
#   return(Q)
# }
# 
# 
# 
# 
# rmvnorm_canonical_Matrix <- function(n = 1, b, Q, Imult = 0){
#   p = ncol(Q)
#   if("sparseMatrix" %in% is(Q) & !("diagonalMatrix" %in% is(Q))){
#     A = Matrix::expand(Matrix::Cholesky(Q, Imult = Imult))
#     w = Matrix::solve(A$L, A$P%*%b)
#     mu = Matrix::solve(Matrix::t(A$L), w)
#     z = rnorm(p*n)
#     dim(z) = c(p,n)
#     v = solve(Matrix::t(A$L), z)
#     y = mu + v #
#     y = Matrix::crossprod(A$P,y) ## P' %*% y
#     out = as.matrix(Matrix::t(y))
#   }else if("diagonalMatrix" %in% is(Q)){
#     variance = as.numeric(1/(Matrix::diag(Q) + Imult)) # length p
#     mu = as.numeric(variance*b) # length p
#     #browser()
#     y = rnorm(p*n, mu, sqrt(variance))
#     dim(y) = c(p, n)
#     out = as.matrix(t(y))
#   }else{
#     if(Imult != 0) diag(Q) = diag(Q) + Imult
#     return(spam::rmvnorm.canonical(n, b, Q))
#     # or alternatively....
#     # R = chol(Q)
#     # w = base::backsolve(R, b, transpose = T)
#     # mu = base::backsolve(R, w)
#     # z = rnorm(p*n)
#     # dim(z) = c(p,n)
#     # v = base::backsolve(R, z)
#     # y = mu + v
#   }
#   return(out)
# }
# 
