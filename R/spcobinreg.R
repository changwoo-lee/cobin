#' spatial cobin regression model
#' 
#' Fit Bayesian spatial cobin regression model under canonical link (cobit link) with Markov chain Monte Carlo (MCMC).
#' \deqn{
#'  y(s_{i}) \mid x(s_{i}), u(s_i) \stackrel{ind}{\sim} cobin(x(s_{i})^T\beta + u(s_i), \lambda^{-1}), \quad u(\cdot)\sim GP
#' }
#' for \eqn{i=1,\dots,n}. See \link[cobin]{dcobin} for details on cobin distribution. It currently only supports mean zero GP with exponential covariance
#' \deqn{
#'   cov(u(s_i), u(s_j)) = \sigma_u^2\exp(-\phi_u d(s_i,s_j))
#' }
#' where \eqn{\phi_u} corresponds to inverse range parameter. 
#' 
#' The prior setting can be controlled with "priors" argument. 
#' Prior for regression coefficients are independent normal or t prior centered at 0.
#' "priors" is a named list of:  
#' \itemize{
#' \item beta_intercept_scale, Default 100, the scale of the intercept prior
#' \item beta_scale, Default 100, the scale of nonintercept fixed-effect coefficients
#' \item beta_df, Default Inf, degree of freedom of t prior. If `beta_df=Inf`, it corresponds to normal prior
#' \item lambda_grid, Default 1:70, candidate for lambda (integer)
#' \item lambda_logprior, Default \eqn{p(\lambda)\propto \lambda \Gamma(\lambda+1)/\Gamma(\lambda+5)}, log-prior of lambda. Default choice arises from beta negative binomial distribution; \eqn{(\lambda-1)\mid \psi \sim negbin(2,\psi), \psi\sim Beta(2,2)}. 
#' \item logprior_sigma.sq, Default half-Cauchy on the sd(u) \eqn{=\sigma_u}, log prior of var(u)\eqn{=\sigma_u^2}
#' \item phi_lb, lower bound of uniform prior of \eqn{\phi_u} (inverse range parameter of spatial random effect). Can be same as phi_ub
#' \item phi_ub, lower bound of uniform prior of \eqn{\phi_u} (inverse range parameter of spatial random effect). Can be same as phi_lb
#' }
#'
#' @param formula an object of class "\link[stats]{formula}" or a two-sided linear formula object describing both the fixed-effects and random-effects part of the model; see "\link[lme4]{lmer}" 
#' @param data data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param link character, link function (default "cobit"). Only supports canonical link function "cobit" that is compatible with Kolmogorov-Gamma augmentation. 
#' @param coords a n x 2 matrix of Euclidean coordinates
#' @param NNGP logical, if TRUE, use NNGP prior for the spatial random effects; see \link[spNNGP]{spNNGP}
#' @param contrasts an optional list. See the contrasts.arg of \link[stats]{model.matrix.default}.
#' @param priors a list of prior hyperparameters. See Details
#' @param nngp.control a list of control parameters for NNGP prior (only when NNGP = TRUE). This should be a named list of n.neighbors and ord, with default of 15 and first coordiate-based ordering. See \link[spNNGP]{spNNGP} for details.
#' @param nburn number of burn-in MCMC iterations.
#' @param nsave number of posterior samples. Total MCMC iteration is nburn + nsave*nthin
#' @param nthin thin-in rate. Total MCMC iteration is nburn + nsave*nthin
#' @import spNNGP
#' @import lme4
#' @import fields
#' @import coda
#'
#' @returns Returns list of 
#' \item{post_save}{a matrix of posterior samples (coda::mcmc) with nsave rows}
#' \item{post_u_save}{a matrix of posterior samples (coda::mcmc) of random effects, with nsave rows}
#' \item{loglik_save}{a nsave x n matrix of pointwise log-likelihood values, can be used for WAIC calculation.}
#' \item{priors}{list of hyperprior information}
#' \item{nsave}{number of MCMC samples}
#' \item{t_mcmc}{wall-clock time for running MCMC}
#' \item{t_premcmc}{wall-clock time for preprocessing before MCMC}
#' \item{y}{response vector}
#' \item{X}{fixed effect design matrix}
#' \item{coords}{a n x 2 matrix of Euclidean coordinates}
#' if NNGP = TRUE, also returns
#' \item{nngp.control}{a list of control parameters for NNGP prior}
#' \item{spNNGPfit}{an "NNGP" class with empty samples, placeholder for prediction}
#' 
#' 
#' @export
#'
#' @examples
#' \dontrun{
#'  # please see https://github.com/changwoo-lee/cobin-reproduce/blob/main/Sec6_mmicasestudy/results_main_n949/run_n949.R
#' }
spcobinreg <- function(formula, data, link = "cobit",
                       coords, NNGP = FALSE, contrasts = NULL,
                     priors = list(beta_intercept_scale = 10,
                                   beta_scale = 2.5, beta_df = Inf),
                     nngp.control = list(n.neighbors = 15, ord = order(coords[, 1])),
                     nburn = 1000, nsave = 1000, nthin = 1){
  if(length(lme4::findbars(formula)) > 0) stop("only supports spatial random intercept model")
  if(link != "cobit") stop("only supports cobit link")

  # lm object for design matrix constructions
  # method = "model.frame": model fit is not performed
  temp_lm = lm(formula, data, contrasts = contrasts, method = "model.frame")
  y = model.response(temp_lm)
  X = model.matrix(temp_lm, data = data)

  
  distmat = fields::rdist(coords)
  


  if(any(y < 0) || any(y > 1)) stop("y must be between 0 and 1")

  # default prior settings
  if(is.null(priors$beta_df)){
    priors$beta_df = Inf
  }
  if(is.null(priors$beta_intercept_scale)){
    priors$beta_intercept_scale = 10
  }
  if(is.null(priors$beta_scale)){
    priors$beta_scale = 2.5
  }
  if(is.null(priors$lambda_grid)){
    priors$lambda_grid = 1:70; 
  }
  lambda_grid = priors$lambda_grid;
  if(is.null(priors$lambda_prior)){
    priors$lambda_logprior = log(lambda_grid) + lfactorial(lambda_grid) - lfactorial(lambda_grid+4) + log(36)
  }
  if(any(y == 0) || any(y == 1)){
    if(priors$lambdagrid != 1) stop("y must be strictly between 0 and 1, unless lambda is fixed as 1")
  }

  if(is.null(priors$logprior_sigma.sq)){
    logprior_sigma.sq = function(x) -log(sqrt(x)*(1+x)) - log(pi) # corresponding to half-cauchy on stdev
    priors$logprior_sigma.sq = logprior_sigma.sq
  }
  # phi : inverse range. 3/phi corresponds to effective range
  # default setting for uniform prior thresholds are between 1/4 and 3/4 of maximal dist

  if(is.null(priors$phi_lb)){ #uniform prior, lower bound
    phi_lb = 3/(max(distmat)*3/4); priors$phi_lb = phi_lb
  }
  if(is.null(priors$phi_ub)){
    phi_ub = 3/(max(distmat)/4); priors$phi_ub = phi_ub
  }
  if(priors$phi_lb == priors$phi_ub){
    print("inverse range parameter fixed")
    priors$phi_fixed = TRUE
  } else{
    priors$phi_fixed = FALSE
  }




  if(NNGP){
    # these values are just dummies (not actually used)
    nngpstarting <- list("phi"=mean(c(priors$phi_lb, priors$phi_ub)), "sigma.sq"=5)
    nngptuning <- list("phi"=mean(c(priors$phi_lb, priors$phi_ub)), "sigma.sq"=0.5)
    nngppriors <- list("phi.Unif"=c(priors$phi_lb-1e-14, priors$phi_ub), "sigma.sq.IG"=c(2, 5))

    ybinary = as.numeric(y > 0.5)
    print("running prelim spNNGP")
    spNNGPfit = spNNGP(formula, data = data, coords = coords,
                       method = "latent", family = "binomial", # this is not important
                       n.neighbors = nngp.control$n.neighbors,
                       starting=nngpstarting, tuning=nngptuning,
                       priors=nngppriors, cov.model="exponential",
                       ord = nngp.control$ord,
                       n.samples=0, verbose = F,
                       return.neighbor.info = TRUE)
    print("completed prelim spNNGP")

    # initial phi
    #phi = mean(c(phi_lb, phi_ub))
    #sigma.sq = 1
    ord = spNNGPfit$neighbor.info$ord
    Nlist = spNNGPfit$neighbor.info$n.indx

    #Sigma = sigma.sq * exp(-distmat * phi) # exponential kernel
    #Q = build_Q_exponential(distmat,
    #                        sigma.sq, phi,
    #                        Nlist = spNNGPfit$neighbor.info$n.indx,
    #                        ord = spNNGPfit$neighbor.info$ord)

  }



  if(!NNGP){
    out = fit_cobin_spatial(y = y, X = X, coords = coords, distmat = distmat,
                            priors = priors,
                             nburn = nburn, nsave = nsave, nthin = nthin)
    out$coords = coords
  }else{
    
    out = fit_cobin_spatial_NNGP(y = y, X = X, coords = coords, distmat = distmat,
                                 priors = priors, ord = ord, Nlist = Nlist,
                                 nburn = nburn, nsave = nsave, nthin = nthin)
    out$coords = coords
    out$nngp.control = nngp.control
    out$spNNGPfit = spNNGPfit # for prediction in the future
  }
  return(out)
}


