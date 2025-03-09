
fit_micobin_spatial_NNGP <- function(y, X, coords, distmat, priors, ord = ord, Nlist = Nlist,
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
  psi_ab = priors$psi_ab
  logprior_sigma.sq = priors$logprior_sigma.sq # function
  phi_lb = priors$phi_lb
  phi_ub = priors$phi_ub
  phi_fixed = priors$phi_fixed
  beta_s = c(beta_intercept_scale, rep(beta_scale,p-1))

  # Initialize
  beta = rep(0,p)
  betavar = beta_s^2
  lambda = rep(1,n)
  psi = 0.5
  gamma = rep(2.5,p)
  kappa = rep(1,n)
  sigma.sq = 1
  phi = mean(c(phi_lb, phi_ub))
  Q = build_Q_exponential(distmat,
                          phi,
                          Nlist = Nlist,
                          ord = ord)

  Q = as(as(Q, "generalMatrix"), "CsparseMatrix")
  print(paste0(round(100-(Matrix::nnzero(Q)-n)/(n^2-n)*100,2), "% of off-diagonal entries of precision matrix is zero"))
  Q = spam::as.spam.dgCMatrix(Q)# spam object
  #Qinvlogdet = -determinant(Q, logarithm = TRUE)$modulus
  Q_spamstruct = spam::chol(Q)

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

    # Step 1-2: sample kappa
    kappa = rkgcpp(n, as.numeric(lambda), as.numeric(linpred))


    ZtKappaX = X*kappa#(t(Z*kappa)%*%X)
    XtKappaX = (t(X*kappa)%*%X)

    # Step 2: sample beta, marginalizing out u
    #nnmat_inv = solve( + Matrix::Diagonal(n, kappa))

    #Vuinv_chol = chol(Q + diag(kappa, nrow = n, ncol = n))
    Vuinv_chol = chol(1/sigma.sq*Q + spam::diag.spam(kappa, n, n))# spam.chol.NgPeython
    #nnmat_inv = chol2inv(Vuinv_chol)
    XtSigma_invX = XtKappaX - crossprod(forwardsolve(Vuinv_chol, ZtKappaX))
    XtSigma_invY = crossprod(X, (y - 0.5)*lambda) - t(forwardsolve(Vuinv_chol, ZtKappaX))%*%forwardsolve(Vuinv_chol, (y - 0.5)*lambda)

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

    if(imcmc < start_adapt){
      proposal = sigma.sq_trans + rnorm(1, 0, sqrt(C0))
    }else{
      proposal = sigma.sq_trans + rnorm(1,0, sqrt(Ct))
    }
    sigma.sq_trans_star = proposal[1]; sigma.sq_star = exp(sigma.sq_trans_star)

    linpred_proxy = (lambda*(y - 0.5) - Xbeta*kappa)/kappa


    # cholSig = chol(diag(1/kappa, nrow = n, ncol = n) + solve(1/sigma.sq*Q))
    # logdetSig = 2*sum(log(diag(cholSig)))
    # logweights = -0.5*sum(backsolve(cholSig, linpred_proxy, transpose = T)^2) - 0.5*logdetSig - n/2*log(2*pi)
    #
    # cholSig = chol(diag(1/kappa, nrow = n, ncol = n) + solve(1/sigma.sq_star*Q))
    # logdetSig = 2*sum(log(diag(cholSig)))
    # logweights_star = -0.5*sum(backsolve(cholSig, linpred_proxy, transpose = T)^2) - 0.5*logdetSig - n/2*log(2*pi)

    QplusKappachol = chol(1/sigma.sq*Q + spam::diag.spam(kappa, n, n))
    quadform = crossprod(sqrt(kappa)*linpred_proxy) - crossprod(forwardsolve(QplusKappachol, kappa*linpred_proxy))
    logweights = -0.5*quadform - determinant(QplusKappachol, logarithm = TRUE)$modulus - 0.5*n*log(sigma.sq) # common factor omitted

    QplusKappachol = chol(1/sigma.sq_star*Q + spam::diag.spam(kappa, n, n))
    quadform = crossprod(sqrt(kappa)*linpred_proxy) - crossprod(forwardsolve(QplusKappachol, kappa*linpred_proxy))
    logweights_star = -0.5*quadform - determinant(QplusKappachol, logarithm = TRUE)$modulus - 0.5*n*log(sigma.sq_star)


    acc_ratio = log(sigma.sq_star) - log(sigma.sq) + # log transformation
      logprior_sigma.sq(sigma.sq_star) - logprior_sigma.sq(sigma.sq) + # uniform prior on rho
      logweights_star - logweights


    if(log(runif(1)) < acc_ratio){
      sigma.sq = sigma.sq_star; sigma.sq_trans = log(sigma.sq)
      Vuinv_chol = chol(1/sigma.sq*Q + spam::diag.spam(kappa, n, n))

      acc_save[imcmc] = 1
      logweights = logweights_star
    }

    # mut and Ct are recursively updated

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

    Vuinv = 1/sigma.sq*Q + spam::diag.spam(kappa, n, n) #diag(kappa, nrow = n, ncol = n) # same as t(D)%*%diag(omega)%*%D
    #b_u = (lambda *(y-0.5) - kappa*Xbeta)
    #u = as.numeric(spam::rmvnorm.canonical(1, b_u, Vuinv))

    b_u = (lambda *(y-0.5) - kappa*Xbeta)
    # Rstruct reduce time a lot
    u = as.numeric(spam::rmvnorm.canonical(1, b_u, Vuinv, Rstruct = Q_spamstruct))
    #
    # # Alg 2.5 of Rue book
    # mutemp = forwardsolve(Vuinv_chol, forwardsolve(Vuinv_chol, b_u, transpose = T))
    # u = mutemp +forwardsolve(Vuinv_chol, rnorm(n))

    #sample psi
    psi = rbeta(1, psi_ab[1] + 2*n, psi_ab[2] - n + sum(lambda))


    # save
    if((imcmc > nburn)&&((imcmc-nburn)%%nthin==0)){
      beta_save[isave,] = beta
      u_save[isave,] = u
      sigma.sq_save[isave,] = sigma.sq
      phi_save[isave,] = phi
      psi_save[isave,] = psi
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
