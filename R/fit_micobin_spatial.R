
fit_micobin_spatial <- function(y, X, coords, distmat, priors,
                              nburn = 100, nsave = 1000, nthin = 1, verbose =TRUE){

  t_start = Sys.time()
  #############################################
  n = nrow(X)
  p = ncol(X)

  # set hyperparameters
  beta_df = priors$beta_df
  beta_intercept_scale = priors$beta_intercept_scale
  beta_scale = priors$beta_scale
  lambda_max = priors$lambda_max

  logprior_sigma.sq = priors$logprior_sigma.sq # function
  phi_lb = priors$phi_lb
  phi_ub = priors$phi_ub
  phi_fixed = priors$phi_fixed
  psi_ab = priors$psi_ab

  beta_s = c(beta_intercept_scale, rep(beta_scale,p-1))


  lambda_max = 70
  # Initialize
  beta = rep(0,p)
  betavar = beta_s^2
  lambda = rep(1,n)
  psi = 0.5
  gamma = rep(2.5,p)
  kappa = rep(1,n)
  sigma.sq = 1
  phi = mean(c(phi_lb, phi_ub))
  Sigma = sigma.sq * exp(-distmat * phi) # exponential kernel
  Sigmainv = chol2inv(chol(Sigma))

  u = rep(0, n)
  Xbeta = X%*%beta

  linpred = Xbeta + u
  # Saving objects
  nmcmc = nburn + nsave*nthin
  beta_save = array(0, dim = c(nsave, p))
  u_save = array(0, dim = c(nsave, n))
  psi_save = matrix(0, nsave, 1)
  sigma.sq_save = matrix(0, nsave, 1)
  phi_save = matrix(0, nsave, 1)
  loglik_save = array(0, dim = c(nsave, n))
  acc_save = numeric(nmcmc)
  # adaptive MH tuning (Harrio)

  MH_eps = 0.001
  if(!phi_fixed) MH_s_d = (2.38)^2/2 else MH_s_d = (2.38)^2# denominator 2 corresponds to dimension (sigu2, rho)
  if(!phi_fixed) C0 = MH_s_d*diag(2) else C0 = MH_s_d
  start_adapt = 100 # adapt after 100 iterations

  # Run MCMC
  # pre-calculate h(y, lambda)
  logh_grid = matrix(0, n, lambda_max)
  for(l in 1:lambda_max){
    logh_grid[,l]  = log(l) + dIH(l*y, l, log = T)
  }
  lvec = 1:lambda_max

  # pre-calculate
  Xtym0.5 = crossprod(X, (y - 0.5))
  Ztym0.5 = (y - 0.5)#crossprod(Z, (y - 0.5))
  # initialize
  ZtKappaX = X*kappa#(t(Z*kappa)%*%X)
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
    linpred = Xbeta + u
    temp = (linpred*y - bft(linpred) + log(1-psi)) %*% t(lvec)
    lambda_logprobs = temp + logh_grid - matrix(log(1-psi), n, lambda_max) + matrix(log(lvec), n, lambda_max, byrow = T) + matrix(2*log(psi), n, lambda_max) # last term does not matter but included for log-likelihood calculation
    loglik = matrixStats::rowLogSumExps(lambda_logprobs) # sum over lambda
    lambda_probs = exp(lambda_logprobs - loglik) # logsumexp
    lambda = sample.rowwise(lambda_probs)
    #for(i in 1:n) lambda[i] = sample(1:lambda_max, 1, prob = lambda_probs[i,])
    #if(any(lambda==1)) browser()

    # Step 1-2: sample kappa
    #t_startomega = Sys.time()
    #for(i in 1:n) kappa[i] = rkg.gamma(1, lambda[i], Xbeta[i])
    kappa = rkgcpp(n, as.numeric(lambda), as.numeric(linpred))


    ZtKappaX = X*kappa#(t(Z*kappa)%*%X)
    XtKappaX = (t(X*kappa)%*%X)

    # Step 2: sample beta, marginalizing out u
    #nnmat_inv = solve( + Matrix::Diagonal(n, kappa))
    Vuinv_chol = chol(Sigmainv + diag(kappa, nrow = n, ncol = n)) # this is bottleneck
    #nnmat_inv = solve(Sigmainv + diag(kappa, nrow = n, ncol = n))
    nnmat_inv = chol2inv(Vuinv_chol)
    XtSigma_invX = XtKappaX - t(ZtKappaX)%*%nnmat_inv%*%ZtKappaX
    XtSigma_invY = crossprod(X, (y - 0.5)*lambda) - t(ZtKappaX)%*%nnmat_inv%*%((y - 0.5)*lambda)
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

    # Step 3-1: sample cov kernel parameters, u marginalized out, beta conditioned on
    # transform sigma.sq to real based on exp transform
    # transform rho \in rho_lb, rho_ub to real based on logistic transform
    # on this transformed space, run random walk with bivariate normal proposal

    sigma.sq_trans = log(sigma.sq)
    if(!phi_fixed) phi_trans = glogit(phi, xmin = phi_lb, xmax = phi_ub)

    if(imcmc < start_adapt){
      if(!phi_fixed){
        proposal = c(sigma.sq_trans, phi_trans) + spam::rmvnorm(1, rep(0,2), C0)
      }else{
        proposal = sigma.sq_trans + rnorm(1, 0, sqrt(C0))
      }
    }else{
      if(!phi_fixed){
        proposal = c(sigma.sq_trans, phi_trans) + spam::rmvnorm(1, rep(0,2), Ct)
      }else{
        proposal = sigma.sq_trans + rnorm(1,0, sqrt(Ct))
      }
    }
    sigma.sq_trans_star = proposal[1]; sigma.sq_star = exp(sigma.sq_trans_star)
    if(!phi_fixed){
      phi_trans_star = proposal[2]; phi_star = inv_glogit(phi_trans_star, phi_lb, phi_ub)
    }else{
      phi_star = phi
    }

    Sigma_star = sigma.sq_star * exp(-distmat * phi_star) # exponential kernel
    linpred_proxy = (lambda*(y - 0.5) - Xbeta*kappa)/kappa


    cholSig = chol(diag(1/kappa, nrow = n, ncol = n) + Sigma)
    logdetSig = 2*sum(log(diag(cholSig)))
    logweights = -0.5*sum(backsolve(cholSig, linpred_proxy, transpose = T)^2) - 0.5*logdetSig - n/2*log(2*pi)

    cholSig = chol(diag(1/kappa, nrow = n, ncol = n) + Sigma_star)
    logdetSig = 2*sum(log(diag(cholSig)))
    logweights_star = -0.5*sum(backsolve(cholSig, linpred_proxy, transpose = T)^2) - 0.5*logdetSig - n/2*log(2*pi)

    if(!phi_fixed){
      acc_ratio = log(sigma.sq_star) - log(sigma.sq) + # log transformation
        (log(phi_star-phi_lb) + log(phi_ub - phi_star)) - (log(phi-phi_lb) + log(phi_ub - phi)) + # log(rho_ub - rho_lb) terms or cancelled out # https://www.wolframalpha.com/input?i=d%2Fdx+log%28%28x-a%29%2F%28b-a%29%2F%281-%28x-a%29%2F%28b-a%29%29%29
        logprior_sigma.sq(sigma.sq_star) - logprior_sigma.sq(sigma.sq) + # uniform prior on rho
        logweights_star - logweights
    }else{
      acc_ratio = log(sigma.sq_star) - log(sigma.sq) + # log transformation
        logprior_sigma.sq(sigma.sq_star) - logprior_sigma.sq(sigma.sq) + # uniform prior on rho
        logweights_star - logweights
    }


    if(log(runif(1)) < acc_ratio){
      sigma.sq = sigma.sq_star; sigma.sq_trans = log(sigma.sq)
      phi = phi_star
      if(!phi_fixed) phi_trans = glogit(phi, xmin = phi_lb, xmax = phi_ub)
      Sigma = Sigma_star
      Sigmainv = chol2inv(chol(Sigma))
      Vuinv_chol = chol(Sigmainv + diag(kappa, nrow = n, ncol = n))

      acc_save[imcmc] = 1
      logweights = logweights_star
    }

    # mut and Ct are recursively updated
    if(!phi_fixed){
      if(imcmc == 1){
        mut = c(sigma.sq_trans, phi_trans)
        Ct = MH_s_d*MH_eps*diag(2)
      }else{
        tmpmu = (mut*(imcmc-1)+c(sigma.sq_trans, phi_trans))/imcmc
        # eq (3) of Haario et al. 2001
        Ct = (imcmc-1)*Ct/imcmc+MH_s_d/imcmc*(imcmc*tcrossprod(mut)-
                                                (imcmc+1)*tcrossprod(tmpmu) +
                                                tcrossprod(c(sigma.sq_trans, phi_trans))+
                                                MH_eps*diag(2))
        mut = tmpmu
      }
    }else{
      if(imcmc == 1){
        mut = sigma.sq_trans
        Ct = MH_s_d*MH_eps
      }else{
        tmpmu = (mut*(imcmc-1)+sigma.sq_trans)/imcmc
        # eq (3) of Haario et al. 2001
        Ct = (imcmc-1)*Ct/imcmc+MH_s_d/imcmc*(imcmc*tcrossprod(mut)-
                                                (imcmc+1)*tcrossprod(tmpmu) +
                                                tcrossprod(sigma.sq_trans)+
                                                MH_eps)
        mut = tmpmu
      }
    }





    #Vuinv = Sigmainv + diag(kappa, nrow = n, ncol = n) # same as t(D)%*%diag(omega)%*%D
    #u = as.numeric(spam::rmvnorm.canonical(1, (lambda *(y-0.5) - kappa*Xbeta), Vuinv))
    b_u = (lambda *(y-0.5) - kappa*Xbeta)
    # Alg 2.5 of Rue book
    mutemp = backsolve(Vuinv_chol, backsolve(Vuinv_chol, b_u, transpose = T))
    u = mutemp + backsolve(Vuinv_chol, rnorm(n))

    # step 3: sample psi
    psi = rbeta(1, psi_ab[1] + 2*n, psi_ab[2] - n + sum(lambda))


    # save
    if((imcmc > nburn)&&((imcmc-nburn)%%nthin==0)){
      beta_save[isave,] = beta
      u_save[isave,] = u
      sigma.sq_save[isave,] = sigma.sq
      phi_save[isave,] = phi
      psi_save[isave,]  = psi
      loglik_save[isave,] = as.numeric(loglik)
      isave = isave + 1
    }
  }
  t_end = Sys.time()
  t_mcmc = difftime(t_end, t_start, units = "secs")


  out = list()
  colnames(beta_save) = colnames(X)
  colnames(u_save) = 1:n
  colnames(sigma.sq_save) = "sigma.sq"
  colnames(phi_save) = "phi"
  colnames(psi_save) = "psi"
  out$post_save = coda::mcmc(cbind(beta_save, sigma.sq_save, phi_save, psi_save))
  out$post_u_save = coda::mcmc(u_save)

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
