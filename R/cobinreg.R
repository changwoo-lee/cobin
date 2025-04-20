#' cobin generalized linear (mixed) models
#' 
#' Fit Bayesian cobin regression model under canonical link (cobit link) with Markov chain Monte Carlo (MCMC). 
#' It supports both fixed-effect only model 
#' \deqn{
#'  y_i \mid x_i \stackrel{ind}{\sim} cobin(x_i^T\beta, \lambda^{-1}),
#' }
#' for \eqn{i=1,\dots,n}, and random effect model (v 1.0.x only supports random intercept), e.g. random intercept model
#' \deqn{
#'  y_{ij} \mid x_{ij}, u_i \stackrel{ind}{\sim} cobin(x_{ij}^T\beta + u_i, \lambda^{-1}), \quad u_i\stackrel{iid}{\sim} N(0, \sigma_u^2)
#' }
#' for \eqn{i=1,\dots,n} (group), and \eqn{j=1,\dots,n_i} (observation within group).
#' 
#' The prior setting can be controlled with "priors" argument. 
#' Prior for regression coefficients are independent t or normal prior centered at 0.
#' "priors" is a named list of:  
#' \item{beta_intercept_scale}{Default 100, the scale of the intercept prior}
#' \item{beta_scale}{Default 100, the scale of nonintercept fixed-effect coefficients}
#' \item{beta_df}{Default Inf, degree of freedom of t prior. If Inf, it corresponds to normal prior}
#' \item{lambda_grid}{Default 1:70, candidate for lambda (integer)}
#' \item{lambda_logprior}{Default \eqn{p(\lambda)\propto }}
#' \item{y}{response vector}
#' \item{X}{fixed effect design matrix}
#' if random effect model, also returns
#' \item{post_u_save}{a matrix of posterior samples (coda::mcmc) of random effects}
#' \item{Z}{random effect design matrix}
#' 
#' By default, prior is set as \eqn{\beta \sim N(0, 10^2)} for regression coefficients, \eqn{p(\lambda) \propto }.eqn{\sigma_u^2 \sim IG(1, 1)}
#'
#' @param formula an object of class "\link[stats]{formula}" or a two-sided linear formula object describing both the fixed-effects and random-effects part of the model; see "\link[lme4]{lmer}" 
#' @param data data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param link character, link function (default "cobit"). Only supports canonical link function "cobit" that is compatible with Kolmogorov-Gamma augmentation. 
#' @param contrasts an optional list. See the contrasts.arg of \link[stats]{model.matrix.default}.
#' @param priors a list of prior parameters. See Details
#' @param nburn number of burn-in MCMC iterations.
#' @param nsave number of posterior samples. Total MCMC iteration is nburn + nsave*nthin
#' @param nthin thin-in rate. Total MCMC iteration is nburn + nsave*nthin
#' @import lme4
#' @import Matrix
#' @import coda
#'
#' @returns Returns list of 
#' \item{post_save}{a matrix of posterior samples (coda::mcmc) with nsave rows}
#' \item{priors}{list of hyperprior information}
#' \item{loglik_save}{a matrix of pointwise log-likelihood values with nsave rows}
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
#' 
#' data("GasolineYield", package = "betareg")
#' 
#' # basic model 
#' out1 = cobinreg(yield ~ temp, data = GasolineYield, 
#'                nsave = 2000, link = "cobit")
#' summary(out1$post_save)
#' plot(out1$post_save)
#' 
#' # random intercept model
#' out2 = cobinreg(yield ~ temp + (1 | batch), data = GasolineYield, 
#'                nsave = 2000, link = "cobit")
#' summary(out2$post_save)
#' plot(out2$post_save)
#' 
cobinreg <- function(formula, data, link = "cobit", contrasts = NULL,
                     priors = list(beta_intercept_scale = 100,
                                   beta_scale = 100, beta_df = Inf),
                     nburn = 1000, nsave = 1000, nthin = 1, MH = F, lambda_fixed = NULL){
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
    priors$beta_intercept_scale = 10
  }
  if(is.null(priors$beta_scale)){
    priors$beta_scale = 2.5
  }
  if(is.null(priors$lambda_grid)){
    priors$lambda_grid = 1:70; 
  }
  lambda_grid = priors$lambda_grid;
  if(is.null(priors$lambda_logprior)){
    priors$lambda_logprior = log(lambda_grid) + lfactorial(lambda_grid) - lfactorial(lambda_grid+4) + log(36)
  }
  if(any(y == 0) || any(y == 1)){
    if(priors$lambdagrid != 1) stop("y must be strictly between 0 and 1, unless lambda is fixed as 1")
  }
  if(isZ){
    if(is.null(priors$a_u)) priors$a_u = 1
    if(is.null(priors$b_u)) priors$b_u = 1
  }

  if(!isZ && MH){
    return(fit_cobin_fixedeffect_mh(y = y, X = X, priors = priors,
                             nburn = nburn, nsave = nsave, nthin = nthin))
  }
  if(!is.null(lambda_fixed)){
    if(isZ) stop("lambda_fixed is not supported for mixed effect models")
    return(fit_cobin_fixedeffect_lambdafixed(y = y, X = X, priors = priors,
                             nburn = nburn, nsave = nsave, nthin = nthin, lambda_fixed = lambda_fixed))
  }
  if(!isZ){
    out = fit_cobin_fixedeffect(y = y, X = X, priors = priors,
                             nburn = nburn, nsave = nsave, nthin = nthin)
  }else{
    out = fit_cobin_mixedeffect(y = y, X = X, Z = Z, priors = priors,
                             nburn = nburn, nsave = nsave, nthin = nthin)
  }
  return(out)
}

# #
# #Example 1
# n = 2000
# p = 2
# X = matrix(rnorm(n*p, 0, sd = 3), n, p)
# 
# X = mvnfast::rmvn(n, rep(0, 2), 9*matrix(c(1, 0.9,
#                                          0.9, 1), 2, 2))
# 
# beta = c(0, 1) 
# y = rcobin(n, X%*%beta, rep(5, n))
# 
# df = data.frame(y = y, X = X)
# out1 = cobinreg(y ~ X, data = df, nburn = 1000, nsave = 10000)
# out2 = cobinreg(y ~ X, data = df, MH = T, nburn = 1000, nsave = 10000, nthin = 1)
# 
# plot(out2$acc_save)
# 
# str(out1$post_save)
# mcmcse::multiESS(out1$post_save[,1:3])/as.numeric(out1$t_mcmc)
# #plot(out1$post_save)
# mcmcse::multiESS(out2$post_save[,1:3])/as.numeric(out2$t_mcmc)
# #plot(out2$post_save)
# 
# coda::effectiveSize(out1$post_save)/as.numeric(out1$t_mcmc)
# coda::effectiveSize(out2$post_save)/as.numeric(out2$t_mcmc)
# 
# 
# 
# 
# 
# 
# n = 1000
# p = 3
# X = mvnfast::rmvn(n, rep(0, p), matrix(c(1, 0.9, 0.5,
#                                          0.9, 1, 0.5,
#                                          0.5, 0.5, 1), p, p))
# beta = c(0, 3, -3) # without intercept
# y = rcobin(n, X%*%beta, rep(10, n))
# 
# df = data.frame(y = y, X = X)
# out1 = cobinreg(y ~ X, data = df, nburn = 1000, nsave = 1000)
# out2 = cobinreg(y ~ X, data = df, MH = T, nburn = 10000, nsave = 1000, nthin = 1)
# coda::effectiveSize(out1$post_save)/as.numeric(out1$t_mcmc)
# coda::effectiveSize(out2$post_save)/as.numeric(out2$t_mcmc)
# 
# plot(out2$post_save)
# 
# 
# mean(out2$acc_save[1001:10000])
# 
# 
# # example 2: random intercept model
# 
# # Define parameters
# n <- 1000
# p <- 3
# X <- matrix(rnorm(n * p), n, p)
# beta <- c(0, 3, -3)  # without intercept
# # Define grouping variable (10 groups with n/10 observations each)
# id <- rep(1:10, each = n/10)
# 
# # Simulate random intercepts for each group
# sigma_re <- 1  # standard deviation for random effects
# u <- rnorm(10, mean = 0, sd = sigma_re)
# u = u - mean(u)
# 
# y = rcobin(n, X%*%beta + u[id], rep(10, n))
# # Combine into a data frame
# df = data.frame(y = y, X=X, id = id)
# 
# out = cobinreg(y ~ X + (1|id), data = df, nburn = 10, nsave = 1000)
# 
#  u
#  (out$post_u_save)
#  summary(out$post_save)
# # str(out$post_save)
# #
# bayesplot::mcmc_trace(out$post_save)
# bayesplot::mcmc_intervals(out$post_u_save)
# u
# bayesplot::mcmc_trace(out$post_save)
# bayesplot::mcmc_trace(out$post_u_save)
# 
# #
# #
# #


