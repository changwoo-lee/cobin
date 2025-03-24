

fit_micobin_fixedeffect_psifixed <- function(y, X, Z, priors, psi_fixed, 
                               nburn = 100, nsave = 1000, nthin = 1,  verbose=TRUE){

  t_start = Sys.time()
  #############################################
  n = nrow(X)
  p = ncol(X)

  # set hyperparameters
  beta_df = priors$beta_df
  beta_intercept_scale = priors$beta_intercept_scale
  beta_scale = priors$beta_scale
  lambda_max = priors$lambda_max
  psi_ab = priors$psi_ab
  beta_s = c(beta_intercept_scale, rep(beta_scale,p-1))

  # Initialize
  beta = rep(0,p)
  betavar = beta_s^2
  lambda = rep(1,n)
  kappa = rep(1,n)
  psi = psi_fixed
  Xbeta = X%*%beta
  # Saving objects
  nmcmc = nburn + nsave*nthin
  beta_save = array(0, dim = c(nsave, p))
  psi_save = matrix(0, nsave, 1)
  loglik_save = array(0, dim = c(nsave, n))

  # Run MCMC
  # pre-calculate h(y, lambda)
  logh_grid = matrix(0, n, lambda_max)
  for(l in 1:lambda_max){
    logh_grid[,l]  = log(l) + dIH(l*y, l, log = T)
  }
  lvec = 1:lambda_max


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
    linpred = as.numeric(Xbeta)
    temp = (linpred*y - bft(linpred) + log(1-psi)) %*% t(lvec)
    lambda_logprobs = temp + logh_grid - matrix(log(1-psi), n, lambda_max) + matrix(log(lvec), n, lambda_max, byrow = T) + matrix(2*log(psi), n, lambda_max) # last term does not matter but included for log-likelihood calculation
    loglik = matrixStats::rowLogSumExps(lambda_logprobs) # sum over lambda
    lambda_probs = exp(lambda_logprobs - loglik) # logsumexp
    lambda = sample.rowwise(lambda_probs)

    # Step 1-2: sample kappa
    kappa = rkgcpp(n, as.numeric(lambda), as.numeric(Xbeta))

    # Step 2: sample beta
    Vbetainv = diag(1/betavar, ncol = p) + crossprod(sqrt(kappa)*X) # same as t(D)%*%diag(omega)%*%D
    #betatilde = solve(Vbetainv, lambda * Xt_ym0.5) # Y: binary vector
    beta = as.numeric(spam::rmvnorm.canonical(1, crossprod(X, (y - 0.5)*lambda), Vbetainv))
    Xbeta = X%*%beta
    # update beta variance for mixture prior
    if(!is.infinite(beta_df)){
      gamma = 1/rgamma(p, shape = beta_df/2 + 1/2, rate = beta_s^2*beta_df/2 + beta^2/2)
    }

    # step 3: sample psi
    #psi = rbeta(1, psi_ab[1] + 2*n, psi_ab[2] - n + sum(lambda))


    # save
    if((imcmc > nburn)&&((imcmc-nburn)%%nthin==0)){
      beta_save[isave,] = beta
      psi_save[isave,] = psi
      loglik_save[isave,] = as.numeric(loglik)
      isave = isave + 1
    }
  }
  t_end = Sys.time()
  t_mcmc = difftime(t_end, t_start, units = "secs")


  out = list()
  colnames(beta_save) = colnames(X)
  colnames(psi_save) = "psi"
  out$post_save = coda::mcmc(cbind(beta_save, psi_save))

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
