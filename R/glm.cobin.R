#' Title
#'
#' @param formula
#' @param data
#' @param weights
#' @param subset
#' @param na.action
#' @param start
#' @param etastart
#' @param mustart
#' @param control
#' @param method
#' @param model
#' @param x
#' @param y
#' @param contrasts
#' @param ...
#' @param lambda_list
#' @param link
#' @param verbose
#'
#' @returns
#' @export
#'
#' @examples
glm.cobin <- function (formula, data, weights, subset, na.action, start = NULL,
                       etastart, mustart, control = glm.control(...), method = "glm.fit",
                       model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, ...,
                       lambda_list = 1:70, link = "cobit", verbose = TRUE)
{

  link <- substitute(link)
  fam0 <- do.call("cobinfamily",list(lambda = 1, link = link)) # initial lambda

  mf <- Call <- match.call()
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(mf)
  Terms <- attr(mf, "terms")
  if (method == "model.frame")
    return(mf)
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "numeric")
  X <- if (!is.empty.model(Terms))
    model.matrix(Terms, mf, contrasts)
  else matrix(, NROW(Y), 0)
  w <- model.weights(mf)
  if (!length(w))
    w <- rep(1, nrow(mf))
  else if (any(w < 0))
    stop("negative weights not allowed")
  offset <- model.offset(mf)
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  n <- length(Y)
  if (!missing(method)) {
    if (!exists(method, mode = "function"))
      stop(gettextf("unimplemented method: %s", sQuote(method)),
           domain = NA)
    glm.fitter <- get(method)
  }
  else {
    method <- "glm.fit"
    glm.fitter <- stats::glm.fit
  }
  fit <- glm.fitter(x = X, y = Y, weights = w, start = start,
                    etastart = etastart, mustart = mustart, offset = offset,
                    family = fam0, control = list(maxit = control$maxit,
                                                  epsilon = control$epsilon, trace = control$trace >
                                                    1), intercept = attr(Terms, "intercept") > 0)
  if(!fit$converge) stop("algorithm (stats::glm.fit) not converged")
  #####################
  betahat = fit$coefficients
  # find lambda that minimizes logliklihood
  ll_save = numeric(length(lambda_list))
  for (l in 1:length(lambda_list)) { # can be parallelized
    ll_save[l] <- sum(dcobin(Y, bftprimeinv(fit$family$linkinv(X%*%betahat)), lambda_list[l], log = T))
  }
  # find lambda
  lambdahat <- lambda_list[which.max(ll_save)]
  if(verbose) print(paste0("Among lambda_list (default 1,...,70), MLE of lambda (inverse dispersion parameter) = ",lambdahat))
  ####################
  fam <- do.call("cobinfamily",list(lambda = lambdahat, link = link))
  ####################
  if (any(Y == 0) || any(Y==1)) warning("Response has bounary values (0 or 1), results should be interpreted as quasi-MLE unless lambda = 1")
  
  fit <- glm.fitter(x = X, y = Y, weights = w, start = start,
                    etastart = etastart, mustart = mustart, offset = offset,
                    family = fam, control = list(maxit = control$maxit,
                                                 epsilon = control$epsilon, trace = control$trace >
                                                   1), intercept = attr(Terms, "intercept") > 0)
  if(!fit$converge) stop("algorithm (stats::glm.fit) not converged")
  Call$lambda <- lambdahat
  Call$link <- link
  #fit$call <- Call

  if (model)
    fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if (x)
    fit$x <- X
  if (!y)
    fit$y <- NULL
  out = structure(c(fit, list(call = Call, formula = formula, terms = mt,
                              data = data, offset = offset, control = control, method = method,
                              contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt,
                                                                                      mf))), class = c(fit$class, c("glm", "lm")))

  out$lambda = lambdahat
  out$lambda_list = lambda_list
  out$ll_save = ll_save
  names(out$ll_save) = paste("lambda =",out$lambda_list)
  return(out)
}

# 
# andrew <- read.csv("/Users/clee/Dropbox/research/cobin/cobinreg/demo/andrew.csv", colClasses = c(QUAD = "factor", PATCH = "factor"))
# andrew$PATCHm <- as.factor(as.numeric(as.character(andrew$PATCH)) %% 5 + 1)
# library(plyr)
# andrew2 <- ddply(andrew, ~PATCH, summarise, ALGAE.mean = mean(ALGAE) / 100, treat = TREAT[1])
# andrew2
# 
# 
# andrew2$treat1234 = as.numeric(as.factor(andrew2$treat))
# 
# 
# out = glm.cobin(ALGAE.mean ~ treat1234, data = andrew2, link = "cobit")
# 
# model.matrix(out)
# 
# 
# summary(out)
# plot(out)
# 







