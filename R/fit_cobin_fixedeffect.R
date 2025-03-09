
fit_cobin_fixedeffect <- function(y, X, Z, priors,
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
  lambda_prior = priors$lambda_prior
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

  # Run MCMC
  # pre-calculate h(y, lambda)
  logh_grid = matrix(0, n, length(lambda_grid))
  for(l in 1:length(lambda_grid)){
    logh_grid[,l]  = log(lambda_grid[l]) + dIH(lambda_grid[l]*y, lambda_grid[l], log = T)
  }
  colsum_logh_grid = colSums(logh_grid)
  Xt_ym0.5 = crossprod(X, y-0.5)

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
    lambda_logprobs = colSums(temp) + colsum_logh_grid + log(lambda_prior)
    lambda = sample(lambda_grid, size = 1, prob = exp(lambda_logprobs - matrixStats::logSumExp(lambda_logprobs) ))

    # Step 1-2: sample kappa
    kappa = rkgcpp(n, as.numeric(rep(lambda, n)), as.numeric(Xbeta))

    # Step 2: sample beta
    Vbetainv = diag(1/betavar, ncol = p) + crossprod(sqrt(kappa)*X) # same as t(D)%*%diag(omega)%*%D
    #betatilde = solve(Vbetainv, lambda * Xt_ym0.5) # Y: binary vector
    beta = as.numeric(spam::rmvnorm.canonical(1, lambda * Xt_ym0.5, Vbetainv))
    Xbeta = X%*%beta
    # update beta variance for mixture prior
    if(!is.infinite(beta_df)){
      gamma = 1/rgamma(p, shape = beta_df/2 + 1/2, rate = beta_s^2*beta_df/2 + beta^2/2)
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

  out$priors = priors

  out$t_mcmc = t_mcmc
  out$t_premcmc = t_premcmc
  out$y = y
  out$X = X
  #out$times = c(t1,t2,t3,t4)
  return(out)
}
