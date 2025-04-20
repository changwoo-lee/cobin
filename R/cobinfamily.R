
# cobit: canonical link for cobin regression
# cobit <- structure(list(linkfun = bftprimeinv, # g
#                         linkinv = bftprime, # ginv
#                         mu.eta = bftprimeprime, #ginvprime
#                         d2mu.deta = bftprimeprimeprime,
#                         valideta = function(eta) TRUE,
#                         name = "cobit"), class = "link-glm")

# similar to MASS negative.binomial

#' cobin family class
#'
#' Specifies the information required to fit a cobin generalized linear model with known lambda parameter, using glm().
#'
#' @param lambda The known value of lambda, must be integer
#' @param link The link function to be used. Options are "cobit" (canonical link for cobin regression), "logit", "probit", "cauchit", "cloglog"
#'
#' @return An object of class "family", a list of functions and expressions needed by glm() to fit a cobin generalized linear model.
#' @export
#'
cobinfamily <- function (lambda = stop("'lambda' must be specified"), link = "cobit")
{
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  if (linktemp %in% c("logit", "probit", "cloglog", "cauchit"))
    stats <- make.link(linktemp)
  else if (linktemp == "cobit"){
    stats = structure(list(linkfun = bftprimeinv, # g
                           linkinv = bftprime, # ginv
                           mu.eta = bftprimeprime, #ginvprime
                           d2mu.deta = bftprimeprimeprime,
                           valideta = function(eta) TRUE,
                           name = "cobit"), class = "link-glm")
  }
  # else if (is.character(link)) {
  #   stats <- make.link(link)
  #   linktemp <- link
  # }
  else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name))
        linktemp <- stats$name
    }
    else stop(gettextf("\"%s\" link not available; available links are \"cobit\" (canonical link for cobin regression), \"logit\", \"probit\", \"cauchit\", \"cloglog\"",
                       linktemp))
  }
  .Lambda <- lambda
  env <- new.env(parent = .GlobalEnv)
  assign(".Lambda", lambda, envir = env)
  variance <- function(mu) Vft(mu)/.Lambda
  validmu <- function(mu) all(mu > 0 & mu < 1)
  dev.resids <- function(y, mu, wt) 2 * wt * (y * ( bftprimeinv(y) - bftprimeinv(mu) ) - bft(bftprimeinv(y)) + bft(bftprimeinv(mu)) )
  aic <- function(y, n, mu, wt, dev) {
    -2 * sum(dcobin(y, bftprimeinv(mu), .Lambda, log = TRUE) * wt)
  }
  initialize <- expression({
    if (any(y < 0)) stop("negative values not allowed")
    if (any(y > 1)) stop("values greater than 1 not allowed")
    n <- rep(1, nobs)
    y[y>0.9999] = 0.9999 # This is just for initializaton
    y[y<0.0001] = 0.0001 # This is just for initializaton
    mustart <- y
  })
  simfun <- function(object, nsim) {
    ftd <- fitted(object)
    rcobin(nsim * length(ftd), bftprimeinv(ftd), .Lambda)
  }
  environment(variance) <- environment(validmu) <- environment(dev.resids) <- environment(aic) <- environment(simfun) <- env
  famname <- paste("cobin(", format(lambda), ")", sep = "")
  structure(list(family = famname, link = linktemp, linkfun = stats$linkfun,
                 linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids,
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
                 validmu = validmu, valideta = stats$valideta, simulate = simfun),
            class = "family")
}

