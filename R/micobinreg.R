#' micobin generalized linear (mixed) models
#' 
#' Fit Bayesian micobin regression model under canonical link (cobit link) with Markov chain Monte Carlo (MCMC). 
#' It supports both fixed-effect only model 
#' \deqn{
#'  y_i \mid x_i \stackrel{ind}{\sim} micobin(x_i^T\beta, \psi),
#' }
#' for \eqn{i=1,\dots,n}, and random intercept model (v 1.0.x only supports random intercept),
#' \deqn{
#'  y_{ij} \mid x_{ij}, u_i \stackrel{ind}{\sim} micobin(x_{ij}^T\beta + u_i, \psi), \quad u_i\stackrel{iid}{\sim} N(0, \sigma_u^2)
#' }
#' for \eqn{i=1,\dots,n} (group), and \eqn{j=1,\dots,n_i} (observation within group). See \link[cobin]{dmicobin} for details on micobin distribution.
#' 
#' The prior setting can be controlled with "priors" argument. 
#' Prior for regression coefficients are independent normal or t prior centered at 0.
#' "priors" is a named list of:  
#' \itemize{
#' \item beta_intercept_scale, Default 100, the scale of the intercept prior
#' \item beta_scale, Default 100, the scale of nonintercept fixed-effect coefficients
#' \item beta_df, Default Inf, degree of freedom of t prior. If `beta_df=Inf`, it corresponds to normal prior
#' \item lambda_max, Default 70, upper bound for lambda (integer)
#' \item psi_ab, Default c(2,2), beta shape parameters for \eqn{\psi} (length 2 vector).
#' }
#' if random intercept model, u ~ InvGamma(a_u,b_u) with
#' \itemize{
#' \item a_u, Default 1, first parameter of Inverse Gamma prior of u
#' \item b_u, Default 1, second parameter of Inverse Gamma prior of u
#' }
#'
#' @param formula an object of class "\link[stats]{formula}" or a two-sided linear formula object describing both the fixed-effects and random-effects part of the model; see "\link[lme4]{lmer}" 
#' @param data data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param link character, link function (default "cobit"). Only supports canonical link function "cobit" that is compatible with Kolmogorov-Gamma augmentation. 
#' @param contrasts an optional list. See the contrasts.arg of \link[stats]{model.matrix.default}.
#' @param priors a list of prior hyperparameters. See Details
#' @param nburn number of burn-in MCMC iterations.
#' @param nsave number of posterior samples. Total MCMC iteration is nburn + nsave*nthin
#' @param psi_fixed logical, fixing psi; experimental
#' @param nthin thin-in rate. Total MCMC iteration is nburn + nsave*nthin
#'
#' @import lme4
#' @import Matrix
#' @import coda
#'
#' @returns Returns list of 
#' \item{post_save}{a matrix of posterior samples (coda::mcmc) with nsave rows}
#' \item{loglik_save}{a nsave x n matrix of pointwise log-likelihood values, can be used for WAIC calculation.}
#' \item{priors}{list of hyperprior information}
#' \item{nsave}{number of MCMC samples}
#' \item{t_mcmc}{wall-clock time for running MCMC}
#' \item{t_premcmc}{wall-clock time for preprocessing before MCMC}
#' \item{y}{response vector}
#' \item{X}{fixed effect design matrix}
#' if random effect model, also returns
#' \item{post_u_save}{a matrix of posterior samples (coda::mcmc) of random effects}
#' \item{Z}{random effect design matrix}
#' 
#' 
#' @export
#'
#' @examples
#' \donttest{
#' requireNamespace("betareg", quietly = TRUE)
#' library(betareg)# for dataset example
#' data("GasolineYield", package = "betareg")
#' 
#' # basic model 
#' out1 = micobinreg(yield ~ temp, data = GasolineYield, 
#'                nsave = 2000, link = "cobit")
#' summary(out1$post_save)
#' plot(out1$post_save)
#' 
#' # random intercept model
#' out2 = micobinreg(yield ~ temp + (1 | batch), data = GasolineYield, 
#'                nsave = 2000, link = "cobit")
#' summary(out2$post_save)
#' plot(out2$post_save)
#' }
micobinreg <- function(formula, data, link = "cobit", contrasts = NULL,
                     priors = list(beta_intercept_scale = 100,
                                   beta_scale = 100, beta_df = Inf),
                     nburn = 1000, nsave = 1000, nthin = 1, psi_fixed = NULL){
  if(link != "cobit") stop("only supports cobit link")

  isZ = length(lme4::findbars(formula)) > 0
  if(!isZ){
    # lm object for design matrix constructions
    # method = "model.frame": model fit is not performed
    temp_lm = lm(formula, data, contrasts = contrasts, method = "model.frame")
    y = model.response(temp_lm)
    X = model.matrix(temp_lm, data = data)
  }else{
    # version 1.0.x: only supports random intercept model
    # TODO: general mixed model with non-diagonal covariance of random effect
    if(as.character((lme4::findbars(formula))[[1]])[2]!="1") stop("version 1.0.x only supports random intercept model")
    
    # lmer object for design matrix constructions
    # optimizer = NULL: model fit is not performed
    temp_lmer = lmer(formula, data, contrasts = contrasts,
                     control = lmerControl(optimizer = NULL))
    if(length(lme4::getME(temp_lmer, "flist"))>1) stop("only support one grouping variable")
    y = lme4::getME(temp_lmer, "y") # response
    X = lme4::getME(temp_lmer, "X") # fixed effect design matrix
    Z = lme4::getME(temp_lmer, "Z") # random effect design matrix
    colnames(Z) = paste0(colnames(Z),colnames(ranef(temp_lmer)[[1]]))
    id = as.numeric(lme4::getME(temp_lmer, "flist")[[1]]) # random effect groups (Assuming one rand effect)
  }

  if(any(y < 0) || any(y > 1)) stop("y must be between 0 and 1")

  # default prior settings
  if(is.null(priors$beta_df)){
    priors$beta_df = Inf
  }
  if(is.null(priors$beta_intercept_scale)){
    priors$beta_intercept_scale = 100
  }
  if(is.null(priors$beta_scale)){
    priors$beta_scale = 100
  }
  if(is.null(priors$lambda_max)){
    priors$lambda_max = 70
  }
  if(is.null(priors$psi_ab)){
    priors$psi_ab = c(2,2)
  }

  if(isZ){
    if(is.null(priors$a_u)) priors$a_u = 1
    if(is.null(priors$b_u)) priors$b_u = 1
  }
  if(!is.null(psi_fixed)){
    if(isZ) stop("psi_fixed is not supported for mixed effect models")
    return(fit_micobin_fixedeffect_psifixed(y = y, X = X, priors = priors,
                                             nburn = nburn, nsave = nsave, nthin = nthin, psi_fixed = psi_fixed))
  }
  if(!isZ){
    out = fit_micobin_fixedeffect(y = y, X = X, priors = priors,
                             nburn = nburn, nsave = nsave, nthin = nthin)
  }else{
    out = fit_micobin_mixedeffect(y = y, X = X, Z = Z, priors = priors,
                                          nburn = nburn, nsave = nsave, nthin = nthin)
  }
  return(out)
}


