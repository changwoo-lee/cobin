
#' micobin generalized linear (mixed) models
#'
#'
#'
#'
#'
#' @param formula
#' @param data
#' @param link
#' @param contrasts
#' @param priors
#' @param nburn
#' @param nsave
#' @param nthin
#' @import lme4
#' @import Matrix
#'
#' @returns
#' @export
#'
#' @examples
micobinreg <- function(formula, data, link = "cobit", contrasts = NULL,
                     priors = list(beta_intercept_scale = 10,
                                   beta_scale = 2.5, beta_df = Inf),
                     nburn = 1000, nsave = 1000, nthin = 1){
  if(link != "cobit") stop("only supports cobit link")

  isZ = length(lme4::findbars(formula)) > 0
  if(!isZ){
    # lm object for design matrix constructions
    # method = "model.frame": model fit is not performed
    temp_lm = lm(formula, data, contrasts = contrasts, method = "model.frame")
    y = model.response(temp_lm)
    X = model.matrix(temp_lm, data = data)
  }else{
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
  if(is.null(priors$lambda_max)){
    priors$lambda_max = 70
  }
  if(is.null(priors$psi_ab)){
    priors$psi_ab = c(2,2)
  }

  if(any(y == 0) || any(y == 1)){
    if(priors$lambdagrid != 1) stop("y must be strictly between 0 and 1, unless lambda is fixed as 1")
  }
  if(isZ){
    if(is.null(priors$a_u)) priors$a_u = 1
    if(is.null(priors$b_u)) priors$b_u = 1
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

#
# source("r/micobin.r")
# # Example 1
# n = 200
# p = 3
# X = matrix(rnorm(n*p), n, p)
# beta = c(0, 3, -3) # without intercept
# y = rmicobin(n, X%*%beta, rep(0.5, n))
# plot(X[,1],y, ylim = c(0,1))
# plot(X[,2],y, ylim = c(0,1))
# plot(X[,3],y, ylim = c(0,1))
# df = data.frame(y = y, X = X)
# out = micobinreg(y ~ X, data = df)
# summary(out$post_save)
# 
# #
# # # example 2: random intercept model
# #
# # # Define parameters
# n <- 5000
# p <- 2
# X <- matrix(rnorm(n * p), n, p)
# beta <- c(-1,1)  # without intercept
# # Define grouping variable (10 groups with n/10 observations each)
# id <- rep(1:10, each = n/10)
# 
# # Simulate random intercepts for each group
# sigma_re <- 1  # standard deviation for random effects
# u <- rnorm(10, mean = 0, sd = sigma_re)
# u = u - mean(u)
# 
# y = rmicobin(n, X%*%beta + u[id], rep(0.5, n))
# # Combine into a data frame
# df = data.frame(y = y, X=X, id = id)
# 
# X%*%beta + u[id]
# 
# out = micobinreg(y ~ X + (1|id), data = df, nburn = 10, nsave = 1000)
# 
# u
# colMeans(out$post_u_save)
# # summary(out$post_save)
# # str(out$post_save)
# #
#  bayesplot::mcmc_intervals(out$post_save)
#  bayesplot::mcmc_trace(out$post_save)
#  bayesplot::mcmc_trace(out$post_u_save)
# #
#  bayesplot::mcmc_intervals(out$post_u_save)
# u



