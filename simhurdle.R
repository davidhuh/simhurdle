### Description:
##    Simulate zero-altered response data based on a hurdle model
##
###  Original Author: David Huh
##
###  Dependencies: aster, MASS
##
###  Arguments:           data = data frame
##              formula.fe.bin = formula for the fixed effect equation (binary)
##              formula.fe.cnt = formula for the fixed effect equation (count)
##              formula.re.bin = formula for the random effect equation (binary)
##              formula.re.cnt = formula for the random effect equation (count)
##                 coef.fe.bin = mean estimates of fixed effects (binary)
##                 coef.fe.cnt = mean estimates of fixed effects (count)
##                 vcov.fe.bin = the variance-covariance matrix of the fixed effects (binary)
##                 vcov.fe.cnt = the variance-covariance matrix of the fixed effects (count)
##                 vcov.re.bin = a variance-covariance matrix of the random effects (binary)
##                 vcov.re.cnt = a variance-covariance matrix of the random effects (count)
##                    dist.cnt = zero-truncated count distribution ("xtpois","xtnb")
##                      od.cnt = mean estimate of negative binomial overdispersion parameter
##                    link.cnt = link function ("log")
##
###  Values:  a vector of simulated zero-altered count responses
##
#

simhurdle <- function(data,
                      formula.fe.bin, formula.fe.cnt, formula.re.bin, formula.re.cnt,
                      coef.fe.bin, coef.fe.cnt,
                      vcov.fe.bin, vcov.fe.cnt, vcov.re.bin, vcov.re.cnt,
                      dist.cnt="xtnb",
                      od.cnt, link.bin="logit", link.cnt="log", idvar="id") {
  
  ## Load libraries for random generation from the random effect and outcome distributions
  require(aster, quietly=TRUE)
  require(MASS, quietly=TRUE)
  
  ## validate required arguments
  if (missing(data))
    stop("Missing a dataset w/ covariate values") 
  if (missing(formula.fe.bin) | missing(formula.fe.bin))
    stop("Missing fixed effect formula")
  if (missing(formula.re.bin) | missing(formula.re.bin))
    stop("Missing random effect formula")
  if (missing(coef.fe.bin) | missing(coef.fe.cnt))
    stop("Missing mean estimates of the fixed effects") 
  if (missing(vcov.fe.bin) | missing(vcov.fe.cnt))
    stop("Missing a variance-covariance matrix for the fixed effects")
  if (missing(vcov.re.bin) | missing(vcov.re.cnt))
    stop("Missing a variance-covariance matrix for the random effects")
  
  ## utility to transform linear predictor to units of the outcome
  unlink <- function(mu, link) {
    if (link=="logit") {
      mu.t <- 1/(1+exp(-mu))
    } else if (link=="log") {
      mu.t <- exp(mu)
    } else mu.t <- mu
    return(mu.t)
  }
  
  ## generate design matrices
  
  # fixed effects
  X.bin.mat <- model.matrix(formula.fe.bin, data=data)
  X.cnt.mat <- model.matrix(formula.fe.cnt, data=data)
  
  # random effects
  Z.bin.mat <- model.matrix(formula.re.bin, data=data)
  Z.cnt.mat <- model.matrix(formula.re.cnt, data=data)
  
  ## simulate random effects
  
  # Number of observations
  N.obs <- nrow(data)
  
  # Number of clusters
  id <- eval(parse(text=(paste("data","$",idvar,sep=""))))
  n.id <- length(unique(id))
  
  sim.re.bin <- mvrnorm(n=n.id, mu=rep(0, ncol(vcov.re.bin)), Sigma=vcov.re.bin)
  sim.re.cnt <- mvrnorm(n=n.id, mu=rep(0, ncol(vcov.re.cnt)), Sigma=vcov.re.cnt)
  
  # duplicate random effects for all observations (days) within individual
  sim.re.bin <- cbind(sim.re.bin[id, ])
  sim.re.cnt <- cbind(sim.re.cnt[id, ])
  
  ## calculate fixed and random effects components of linear predictor
  mu.fe.bin <- X.bin.mat %*% coef.fe.bin
  mu.fe.cnt <- X.cnt.mat %*% coef.fe.cnt
  
  mu.re.bin <- rep(NA, N.obs)
  mu.re.cnt <- rep(NA, N.obs)
  
  for (i in 1:N.obs) {
    mu.re.bin[i] <- Z.bin.mat[i,] %*% sim.re.bin[i,]
    mu.re.cnt[i] <- Z.cnt.mat[i,] %*% sim.re.cnt[i,]
  }
    
  mu.bin <- unlink(mu=(mu.fe.bin + mu.re.bin), link=link.bin)
  mu.cnt <- unlink(mu=(mu.fe.cnt + mu.re.cnt), link=link.cnt)
  
  
  ## simulate from the zero-truncated count distribution
  y.sim.cnt <- rep(NA, N.obs)
  
  if (dist.cnt=="xtpois") {  # Zero-truncated Poisson
    for (i in 1:N.obs) {
      y.sim.cnt[i] <- rktp(n=1, k=0 , mu=mu.cnt[i], xpred=1)
    }  
  }
  
  if (dist.cnt=="xtnb") {     # Zero-truncated Negative Binomial
    for (i in 1:N.obs) {
      y.sim.cnt[i] <- rktnb(n=1, size=od.cnt, k=0 , mu=mu.cnt[i], xpred=1)
    }  
  }
  
  ## simulate from binomial distribution
  y.sim.bin <- rep(NA, N.obs)
  for (i in 1:N.obs) {
    y.sim.bin[i] <- rbinom(n=1, size=1, prob=mu.bin[i])
  }
    
  ## combine count and binary simulations
  y.sim <- ifelse(y.sim.bin==0, 0, y.sim.cnt)
  
  ## return vector of responses
  res <- y.sim
  return(res)
}
