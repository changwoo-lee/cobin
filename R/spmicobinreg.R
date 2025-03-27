



#' Title
#'
#' @param formula
#' @param data
#' @param link
#' @param coords
#' @param NNGP
#' @param contrasts
#' @param priors
#' @param nngp.control
#' @param nburn
#' @param nsave
#' @param nthin
#'
#' @returns
#' @export
#'
#' @examples
spmicobinreg <- function(formula, data, link = "cobit",
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
  if(is.null(priors$lambda_max)){
    priors$lambda_max = 70
  }
  if(is.null(priors$psi_ab)){
    priors$psi_ab = c(2,2)
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
    spNNGPfit = spNNGP(ybinary ~ X, data = data, coords = coords,
                       method = "latent", family = "binomial", # this is not important
                       n.neighbors = nngp.control$n.neighbors,
                       starting=nngpstarting, tuning=nngptuning,
                       priors=nngppriors, cov.model="exponential",
                       ord = nngp.control$ord,
                       n.samples=0, verbose = FALSE,
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
    out = fit_micobin_spatial(y = y, X = X, coords = coords, distmat = distmat,
                            priors = priors,
                            nburn = nburn, nsave = nsave, nthin = nthin)
    out$coords = coords
  }else{
    out = fit_micobin_spatial_NNGP(y = y, X = X, coords = coords, distmat = distmat,
                                 priors = priors, ord = ord, Nlist = Nlist,
                                 nburn = nburn, nsave = nsave, nthin = nthin)
    out$coords = coords
    out$nngp.control = nngp.control
    out$spNNGPfit = spNNGPfit # for prediction in the future
  }
  return(out)
}


# 
# n <- 1000
# p <- 1
# X <- matrix(rnorm(n * p), n, p)
# beta <- 1  # without intercept
# 
# library(fields)
# coords = cbind(runif(n), runif(n))
# dmat = rdist(coords)
# Sigma = 1*Matern(dmat, range = 1/0.1, smoothness = 0.5)
# u =as.numeric(mvnfast::rmvn(1, rep(0, n), Sigma))
# y = rmicobin(n, X%*%beta + u, rep(0.5, n))
# 
# df = data.frame(y = y, X=X)
# 
# quilt.plot(coords, u - mean(u))
# 
# myout1 = spmicobinreg(y ~ X, df, coords = coords, NNGP = F, priors = list(phi_lb = 0.1, phi_ub = 0.1),
#                     #nngp.control = list(n.neighbors = 30, ord = order(coords[, 1])),
#                     nburn = 100, nsave = 1000, nthin = 1)
# myout1$t_mcmc
# quilt.plot(coords, colMeans(myout1$post_u_save))
# summary(myout1$post_save)
# plot(myout1$post_save)
# summary(myout1$post_save[,"sigma.sq"])
# 
# library(GpGp)
# myord = GpGp::order_maxmin(coords)
# 
# 
# library(spNNGP)
# myout2 = spmicobinreg(y ~ X, df, coords = coords, NNGP = T, priors = list(phi_lb = 0.1, phi_ub = 0.1),
#                     #nngp.control = list(n.neighbors = 30, ord = myord),
#                     nburn = 100, nsave = 1000, nthin = 1)
# myout2$t_mcmc
# quilt.plot(coords, colMeans(myout2$post_u_save))
# summary(myout2$post_save)
# plot(myout2$post_save)
# quilt.plot(coords, u)
# bayesplot::mcmc_intervals(myout1$post_save)
# bayesplot::mcmc_intervals(myout2$post_save)
# #
# 
# 
# 
# 
