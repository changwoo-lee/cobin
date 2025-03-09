



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
spmicobinreg <- function(formula, data, link = "cobit",
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
  if(is.null(priors$lambda_max)){
    priors$lambda_max = 70
  }
  if(is.null(priors$psi_ab)){
    priors$psi_ab = c(1,1)
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
  }else{
    out = fit_micobin_spatial_NNGP(y = y, X = X, coords = coords, distmat = distmat,
                                 priors = priors, ord = ord, Nlist = Nlist,
                                 nburn = nburn, nsave = nsave, nthin = nthin)
  }
  return(out)
}






