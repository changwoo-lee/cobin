
fit_cobin_fixedeffect_mh <- function(y, X, Z, priors,
                                  nburn = 100, nsave = 1000, nthin = 1,  verbose=TRUE){
  
  t_start = Sys.time()
  #############################################
  n = nrow(X)
  p = ncol(X)
  
  # set hyperparameters
  beta_df = priors$beta_df
  beta_intercept_scale = priors$beta_intercept_scale
  beta_scale = priors$beta_scale
  lambda_grid = priors$lambda_grid
  lambda_logprior = priors$lambda_logprior
  beta_s = c(beta_intercept_scale, rep(beta_scale,p-1))
  
  # Initialize
  beta = rep(0,p)
  betavar = beta_s^2
  lambda = lambda_grid[1]
  kappa = rep(1,n)
  Xbeta = X%*%beta
  # Saving objects
  nmcmc = nburn + nsave*nthin
  beta_save = array(0, dim = c(nsave, p))
  lambda_save = matrix(0, nsave, 1)
  loglik_save = array(0, dim = c(nsave, n))
  acc_save = numeric(nmcmc)
  # Run MCMC
  # pre-calculate h(y, lambda)
  logh_grid = matrix(0, n, length(lambda_grid))
  for(l in 1:length(lambda_grid)){
    logh_grid[,l]  = log(lambda_grid[l]) + dIH(lambda_grid[l]*y, lambda_grid[l], log = T)
  }
  colsum_logh_grid = colSums(logh_grid)
  Xt_ym0.5 = crossprod(X, y-0.5)
  
  MH_eps = 0.001
  MH_s_d = (2.38)^2/p 
  C0 = MH_s_d*diag(p) 
  start_adapt = 100 # adapt after 100 iterations
  
  
  t_end = Sys.time()
  t_premcmc = difftime(t_end, t_start, units = "secs")
  if(verbose){
    pb <- txtProgressBar(style=3)
  }
  isave = 1
  ##### Start of MCMC #####
  t_start = Sys.time()
  for(imcmc in 1:(nburn + nsave*nthin)){
    if(verbose){
      setTxtProgressBar(pb, imcmc/(nburn + nsave*nthin))
    }
    
    
    # Step 1-1: sample lambda
    temp = (Xbeta*y - bft(Xbeta)) %*% t(lambda_grid)
    lambda_logprobs = colSums(temp) + colsum_logh_grid + lambda_logprior
    lambda = sample(lambda_grid, size = 1, prob = exp(lambda_logprobs - matrixStats::logSumExp(lambda_logprobs) ))
    
    # Step 1-2: sample beta
    if(imcmc < start_adapt){
      beta_star = beta + as.numeric(spam::rmvnorm(1, rep(0,p), C0))
    }else{
      beta_star = beta + as.numeric(spam::rmvnorm(1, rep(0,p), Ct))
    }
    Xbeta_star = X%*%beta_star
    acc_ratio = #sum(dnorm(beta_star, rep(0,p), beta_s, log = T) - dnorm(beta, rep(0,p), beta_s, log = T)) + # uniform prior on rho
      lambda*sum(y*(Xbeta_star - Xbeta)) +
      lambda*sum(-bft(Xbeta_star) + bft(Xbeta)) # likelihood
    
    if(log(runif(1)) < acc_ratio){
      beta = beta_star
      Xbeta = Xbeta_star
      acc_save[imcmc] = 1
    }
    
    # mut and Ct are recursively updated
    
    if(imcmc == 1){
      mut = beta
      Ct = MH_s_d*MH_eps
    }else{
      tmpmu = (mut*(imcmc-1)+beta)/imcmc
      # eq (3) of Haario et al. 2001
      Ct = (imcmc-1)*Ct/imcmc+MH_s_d/imcmc*(imcmc*tcrossprod(mut)-
                                              (imcmc+1)*tcrossprod(tmpmu) +
                                              tcrossprod(beta)+
                                              MH_eps*diag(p))
      mut = tmpmu
    }
    
    
    # save
    if((imcmc > nburn)&&((imcmc-nburn)%%nthin==0)){
      beta_save[isave,] = beta
      lambda_save[isave,] = lambda
      loglik_save[isave,] = as.numeric(logh_grid[,which(lambda==lambda_grid)] + lambda*Xbeta*y - lambda*bft(Xbeta))
      isave = isave + 1
    }
  }
  t_end = Sys.time()
  t_mcmc = difftime(t_end, t_start, units = "secs")
  
  
  out = list()
  colnames(beta_save) = colnames(X)
  colnames(lambda_save) = "lambda"
  out$post_save = coda::mcmc(cbind(beta_save, lambda_save))
  
  out$loglik_save = loglik_save
  out$nsave = nsave
  out$acc_save = acc_save
  out$Ct = Ct
  out$priors = priors
  
  out$t_mcmc = t_mcmc
  out$t_premcmc = t_premcmc
  out$y = y
  out$X = X
  #out$times = c(t1,t2,t3,t4)
  return(out)
}
