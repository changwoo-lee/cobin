

fit_micobin_mixedeffect <- function(y, X, Z, priors,
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
  a_u = priors$a_u
  b_u = priors$b_u
  beta_s = c(beta_intercept_scale, rep(beta_scale,p-1))

  # Initialize
  beta = rep(0,p)
  betavar = beta_s^2
  uvar = 1
  lambda = rep(1,n)
  gamma = rep(2.5,p)
  kappa = rep(1,n)
  psi = 0.5
  u = rep(0, ncol(Z))
  Xbeta = X%*%beta
  Zu = Z%*%u
  linpred = Xbeta + Zu
  # Saving objects
  nmcmc = nburn + nsave*nthin
  beta_save = array(0, dim = c(nsave, p))
  u_save = array(0, dim = c(nsave, ncol(Z)))
  psi_save = matrix(0, nsave, 1)
  uvar_save = matrix(0, nsave, 1)
  loglik_save = array(0, dim = c(nsave, n))

  # Run MCMC
  # pre-calculate h(y, lambda)
  logh_grid = matrix(0, n, lambda_max)
  for(l in 1:lambda_max){
    logh_grid[,l]  = log(l) + dIH(l*y, l, log = T)
  }
  lvec = 1:lambda_max

  # initialize
  ZtKappaX = (t(Z*kappa)%*%X)
  XtKappaX = (t(X*kappa)%*%X)

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
    linpred = as.numeric(Xbeta + Zu)
    temp = (linpred*y - bft(linpred) + log(1-psi)) %*% t(lvec)
    lambda_logprobs = temp + logh_grid - matrix(log(1-psi), n, lambda_max) + matrix(log(lvec), n, lambda_max, byrow = T) + matrix(2*log(psi), n, lambda_max) # last term does not matter but included for log-likelihood calculation
    loglik = matrixStats::rowLogSumExps(lambda_logprobs) # sum over lambda
    lambda_probs = exp(lambda_logprobs - loglik) # logsumexp
    lambda = sample.rowwise(lambda_probs)

    # Step 1-2: sample kappa
    kappa = rkgcpp(n, as.numeric(lambda), as.numeric(linpred))

    ZtKappaX = (t(Z*kappa)%*%X)
    XtKappaX = (t(X*kappa)%*%X)

    # Step 3: sample beta, marginalizing out u
    nnmat_inv = solve(Matrix::Diagonal(ncol(Z), 1/uvar) + crossprod(sqrt(kappa)*Z))
    XtSigma_invX = XtKappaX - t(ZtKappaX)%*%nnmat_inv%*%ZtKappaX
   # XtSigma_invY = Xtym0.5*lambda - t(ZtKappaX)%*%nnmat_inv%*%Ztym0.5*lambda
    XtSigma_invY = crossprod(X, (y - 0.5)*lambda) - t(ZtKappaX)%*%nnmat_inv%*%crossprod(Z,(y - 0.5)*lambda)

    if(!is.infinite(beta_df)){ # normal prior
      Q_beta = XtSigma_invX + diag(1/gamma, p)
    }else{ # t prior
      Q_beta = XtSigma_invX + diag(1/beta_s^2, p)
    }
    b_beta = XtSigma_invY # assuming prior mean is zero
    beta = as.numeric(spam::rmvnorm.canonical(1, b_beta, Q_beta))

    Xbeta = X%*%beta
    # update beta variance for mixture prior
    if(!is.infinite(beta_df)){
      gamma = 1/rgamma(p, shape = beta_df/2 + 1/2, rate = beta_s^2*beta_df/2 + beta^2/2)
    }

    # Step 4: sample u,
    Vuinv = Matrix::Diagonal(ncol(Z), 1/uvar) + crossprod(sqrt(kappa)*Z) # same as t(D)%*%diag(omega)%*%D
    u = as.numeric(spam::rmvnorm.canonical(1,  crossprod(Z, lambda *(y-0.5) - kappa*Xbeta), Vuinv))
    Zu = Z%*%u

    # update u variance
    #browser()
    #uvar = 1
    uvar = 1/rgamma(1, shape = a_u + ncol(Z)/2, rate = b_u + sum(u^2)/2)

    # step 3: sample psi
    psi = rbeta(1, psi_ab[1] + 2*n, psi_ab[2] - n + sum(lambda))



    # save
    if((imcmc > nburn)&&((imcmc-nburn)%%nthin==0)){
      beta_save[isave,] = beta
      u_save[isave,] = u
      uvar_save[isave] = uvar
      psi_save[isave,] = psi
      loglik_save[isave,] = as.numeric(loglik)
      isave = isave + 1
    }
  }
  t_end = Sys.time()
  t_mcmc = difftime(t_end, t_start, units = "secs")


  out = list()
  colnames(beta_save) = colnames(X)
  colnames(u_save) = colnames(Z)
  colnames(uvar_save) = "var(u)"
  colnames(psi_save) = "psi"
  out$post_save = coda::mcmc(cbind(beta_save, uvar_save, psi_save))
  out$post_u_save = coda::mcmc(u_save)

  out$loglik_save = loglik_save
  out$nsave = nsave

  out$priors = priors

  out$t_mcmc = t_mcmc
  out$t_premcmc = t_premcmc
  out$y = y
  out$X = X
  out$Z = Z
  #out$times = c(t1,t2,t3,t4)
  return(out)
}
