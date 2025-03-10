
fit_cobin_spatial_NNGP <- function(y, X, coords, distmat, priors, ord = ord, Nlist = Nlist,
                                   nburn = 100, nsave = 1000, nthin = 1, verbose =TRUE){
  
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
  
  logprior_sigma.sq = priors$logprior_sigma.sq # function
  phi_lb = priors$phi_lb
  phi_ub = priors$phi_ub
  phi_fixed = priors$phi_fixed
  beta_s = c(beta_intercept_scale, rep(beta_scale,p-1))
  
  # Initialize
  beta = rep(0,p)
  betavar = beta_s^2
  lambda = lambda_grid[1]
  gamma = rep(2.5,p)
  kappa = rep(1,n)
  sigma.sq = 1
  phi = mean(c(phi_lb, phi_ub))
  Q = build_Q_exponential(distmat,
                           sigma.sq, 
                          phi,
                          Nlist = Nlist,
                          ord = ord)
  logdetQ = Matrix::determinant(Q, logarithm = TRUE)$modulus
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
  lambda_save = matrix(0, nsave, 1)
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
  logh_grid = matrix(0, n, length(lambda_grid))
  for(l in 1:length(lambda_grid)){
    logh_grid[,l]  = log(lambda_grid[l]) + dIH(lambda_grid[l]*y, lambda_grid[l], log = T)
  }
  colsum_logh_grid = colSums(logh_grid)
  
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
    temp = (linpred*y - bft(linpred)) %*% t(lambda_grid)
    lambda_logprobs = colSums(temp) + colsum_logh_grid + lambda_logprior
    lambda = sample(lambda_grid, size = 1, prob = exp(lambda_logprobs - matrixStats::logSumExp(lambda_logprobs) ))
    
    # Step 1-2: sample kappa
    kappa = rkgcpp(n, as.numeric(rep(lambda, n)), as.numeric(linpred))
    #kappa = samplekg::rkg(n, as.numeric(rep(lambda, n)), as.numeric(linpred))
    
    
    ZtKappaX = X*kappa#(t(Z*kappa)%*%X)
    XtKappaX = (t(X*kappa)%*%X)
    
    # Step 2: sample beta, marginalizing out u
    #nnmat_inv = solve( + Matrix::Diagonal(n, kappa))
    #Vuinv_chol = chol(Q + diag(kappa, nrow = n, ncol = n))
    
    #QplusKappachol = chol(Q + spam::diag.spam(kappa, n, n), Rstruct = Q_spamstruct)# spam.chol.NgPeython
    QplusKappachol <- spam::update.spam.chol.NgPeyton(Q_spamstruct, Q + spam::diag.spam(kappa, n, n))
    
    #nnmat_inv = chol2inv(Vuinv_chol)
    XtSigma_invX = XtKappaX - crossprod(spam::forwardsolve(QplusKappachol, ZtKappaX))
    XtSigma_invY = Xtym0.5*lambda - t(spam::forwardsolve(QplusKappachol, ZtKappaX))%*%spam::forwardsolve(QplusKappachol, Ztym0.5*lambda)
    
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
    Q_star= build_Q_exponential(distmat,
                             sigma.sq_star, 
                             phi_star,
                             Nlist = Nlist,
                             ord = ord)
    logdetQ_star = Matrix::determinant(Q_star, logarithm = TRUE)$modulus
    Q_star = as(as(Q_star, "generalMatrix"), "CsparseMatrix")
    Q_star = spam::as.spam.dgCMatrix(Q_star)# spam object
    #browser()
    linpred_proxy = (lambda*(y - 0.5) - Xbeta*kappa)/kappa
    
    
    # cholSig = chol(diag(1/kappa, nrow = n, ncol = n) + solve(1/sigma.sq*Q))
    # logdetSig = 2*sum(log(diag(cholSig)))
    # logweights = -0.5*sum(backsolve(cholSig, linpred_proxy, transpose = T)^2) - 0.5*logdetSig - n/2*log(2*pi)
    #
    # cholSig = chol(diag(1/kappa, nrow = n, ncol = n) + solve(1/sigma.sq_star*Q))
    # logdetSig = 2*sum(log(diag(cholSig)))
    # logweights_star = -0.5*sum(backsolve(cholSig, linpred_proxy, transpose = T)^2) - 0.5*logdetSig - n/2*log(2*pi)
    
    #QplusKappachol = chol(Q + spam::diag.spam(kappa, n, n), Rstruct = Q_spamstruct)
    #QplusKappachol = spam::update.spam.chol.NgPeyton(Q_spamstruct, Q + spam::diag.spam(kappa, n, n))
    
    #quadform = crossprod(sqrt(kappa)*linpred_proxy) - crossprod(forwardsolve(QplusKappachol, kappa*linpred_proxy))
    quadform = - crossprod(spam::forwardsolve(QplusKappachol, kappa*linpred_proxy))
    logweights = -0.5*quadform - spam::determinant(QplusKappachol, logarithm = TRUE)$modulus + 0.5*logdetQ# common factor omitted
    
    #QplusKappachol_star = chol(Q_star + spam::diag.spam(kappa, n, n), Rstruct = Q_spamstruct)
    QplusKappachol_star = spam::update.spam.chol.NgPeyton(Q_spamstruct, Q_star + spam::diag.spam(kappa, n, n))
    #quadform = crossprod(sqrt(kappa)*linpred_proxy) - crossprod(forwardsolve(QplusKappachol_star, kappa*linpred_proxy))
    quadform = - crossprod(spam::forwardsolve(QplusKappachol_star, kappa*linpred_proxy))
    logweights_star = -0.5*quadform - spam::determinant(QplusKappachol_star, logarithm = TRUE)$modulus + 0.5*logdetQ_star
  
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
      Q = Q_star
      logdetQ = logdetQ_star
      QplusKappachol = QplusKappachol_star
      #Vuinv_chol = chol(1/sigma.sq*Q + spam::diag.spam(kappa, n, n))
      acc_save[imcmc] = 1
    }
    
    # mut and Ct are recursively updated
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
    
    
    #Vuinv = Q + spam::diag.spam(kappa, n, n) #diag(kappa, nrow = n, ncol = n) # same as t(D)%*%diag(omega)%*%D
    #b_u = (lambda *(y-0.5) - kappa*Xbeta)
    #u = as.numeric(spam::rmvnorm.canonical(1, b_u, Vuinv))
    
    b_u = (lambda *(y-0.5) - kappa*Xbeta)
    # Rstruct reduce time a lot
    #u = as.numeric(spam::rmvnorm.canonical(1, b_u, Vuinv))
    #u = as.numeric(spam::rmvnorm.canonical(1, b_u, Vuinv, Rstruct = Q_spamstruct))
     
    mu <- as.numeric(solve.spam(QplusKappachol, b_u))
    u = mu + as.numeric(spam::backsolve(QplusKappachol, rnorm(n)))
    # # Alg 2.5 of Rue book
    # mutemp = forwardsolve(Vuinv_chol, forwardsolve(Vuinv_chol, b_u, transpose = T))
    # u = mutemp +forwardsolve(Vuinv_chol, rnorm(n))
    
    
    
    # save
    if((imcmc > nburn)&&((imcmc-nburn)%%nthin==0)){
      beta_save[isave,] = beta
      u_save[isave,] = u
      sigma.sq_save[isave,] = sigma.sq
      phi_save[isave,] = phi
      lambda_save[isave,] = lambda
      loglik_save[isave,] = as.numeric(logh_grid[,which(lambda==lambda_grid)] + lambda*Xbeta*y - lambda*bft(Xbeta))
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
  colnames(lambda_save) = "lambda"
  out$post_save = coda::mcmc(cbind(beta_save, sigma.sq_save, phi_save, lambda_save))
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
