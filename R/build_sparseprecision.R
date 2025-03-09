
build_Q_exponential <- function(distmat, phi,
                                Nlist, ord){
  R = exp(-distmat * phi) # exponential kernel
  n = ncol(R)
  # sparse matrix
  A = Matrix::Matrix(0, n, n);
  D = numeric(n) #D = Matrix(0, n, n);
  D[1] = 1
  # pseudocode 2 of Finley et al.
  #Nlist[[i]] be the set of indices j  <i such that A[i,j] not qe 0
  for(i in 1:(n-1)) {
    temp =  solve(R[Nlist[[i+1]], Nlist[[i+1]]],R[Nlist[[i+1]], i+1])
    A[i+1, Nlist[[i+1]]]= temp
    D[i+1]= R[i+1, i+1] - sum(R[i+1, Nlist[[i+1]]]*temp)
  }
  Dmat = Matrix::Diagonal(n, x = D)
  Q_reordered = Matrix::t((Matrix::Diagonal(n) - A))%*%Matrix::solve(Dmat)%*%(Matrix::Diagonal(n) - A)
  #Linv =  Matrix::Diagonal(n, x = 1/sqrt(D))%*%(Matrix::Diagonal(n) - A)
  Q = Q_reordered[order(ord),order(ord)]
  Q = as(Q, "dgCMatrix")
  return(Q)
}
