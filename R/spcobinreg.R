



#' Title
#'
#' @param formula
#' @param data
#' @param link
#' @param coords
#' @param lonlat
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
spcobinreg <- function(formula, data, link = "cobit",
                       coords, lonlat = FALSE, NNGP = FALSE, contrasts = NULL,
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

  if(lonlat){
    distmat = fields::rdist.earth(coords, miles = FALSE)
  }else{
    distmat = fields::rdist(coords)
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
    priors$lambda_grid = 1:70
  }
  if(is.null(priors$lambda_prior)){
    priors$lambda_prior = rep(1,length(priors$lambda_grid)) #or 0.9^priors$lambda_grid
    priors$lambda_prior = priors$lambda_prior/sum(priors$lambda_prior)
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
    spNNGPfit = spNNGP(ybinary ~ X, data = data, coords = coords,
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
  }else{
    out = fit_cobin_spatial_NNGP(y = y, X = X, coords = coords, distmat = distmat,
                                 priors = priors, ord = ord, Nlist = Nlist,
                                 nburn = nburn, nsave = nsave, nthin = nthin)
  }
  return(out)
}





#
#
#
#
#
#
#
# library(spNNGP)
#
#
#
#
#
# # example: spatial model
#
# # Define parameters
# n <- 1000
# p <- 1
# X <- matrix(rnorm(n * p), n, p)
# beta <- 1  # without intercept
#
# library(fields)
# coords = cbind(runif(n), runif(n))
# dmat = rdist(coords)
# Sigma = 1*Matern(dmat, range = 1/0.1, smoothness = 0.5)
# u =as.numeric(rmvnorm(1, rep(0, n), Sigma))
# quilt.plot(coords, u)
# y = rcobin(n, X%*%beta + u, rep(10, n))
#
# df = data.frame(y = y, X=X)
#
# myout = spcobinreg(y ~ X, df, coords = coords, NNGP = T, priors = list(phi_lb = 0.1, phi_ub = 0.1),
#                    #nngp.control = list(n.neighbors = 30, ord = order(coords[, 1])),
#                    nburn = 100, nsave = 1000)
#
# quilt.plot(coords, u)
#
# quilt.plot(coords, colMeans(myout$post_u_save))
# plot(myout$post_save)
#
# ybinary = as.numeric(y>0.5)
#
# # these are just dummies
# n.samples <- 1000
# starting <- list("phi"=3/0.5, "sigma.sq"=5, "tau.sq"=1)
# tuning <- list("phi"=0.5, "sigma.sq"=0.5, "tau.sq"=0.5)
# priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))
# cov.model <- "exponential"
#
# spNNGPfit = spNNGP(y ~ X, data = df, coords = coords, method = "latent", family = "gaussian",
#                    n.neighbors = nngp.control$n.neighbors,
#                    starting=starting, tuning=tuning,
#                    family = "gaussian", n.neighbors = 15, priors=priors, cov.model=cov.model,
#                    n.samples=n.samples, return.neighbor.info = TRUE)
#
# spNNGPfit$n.neighbors
# spNNGPfit$cov.model
#
# str(spNNGPfit$neighbor.info)
#
#
# ntest = 900
# coords_test = expand.grid(x = seq(0,1, length=30), y = seq(0,1, length=30))
# coords_test = as.matrix(coords_test)
# s.pred <- predict(spNNGPfit, X.0=cbind(1,rep(0,ntest)), coords.0=coords_test,
#                   sub.sample=list(start=501, thin=1),
#                   n.report=1000)
#
# str(s.pred)
#
# quilt.plot(coords_test, rowMeans(s.pred$p.w.0))
