
# inverse cdf of continuous Bernoulli (continuous binomial with lambda = 1)
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

# cdf of continuous Bernoulli (continuous binomial with lambda = 1)
pcb <- function(q, theta){
  theta = array(theta,length(q))
  zeroidx = which(abs(theta) < 1e-12)
  nonzeroidx = which((abs(theta) >= 1e-12) & (theta < 500))
  bigidx = which(theta >= 500)
  out = numeric(length(q))
  out[zeroidx] = q[zeroidx]
  out[nonzeroidx] = expm1(theta[nonzeroidx] * q[nonzeroidx])/expm1(theta[nonzeroidx])
  out[bigidx] = exp(theta[bigidx]*(q[bigidx] - 1))
  return(out)
}

# random variate generation of continuous Bernoulli (continuous binomial with lambda = 1)
rcb <- function(n, theta){
  return(qcb(p = runif(n), theta))
}
