# simulate_exp generates p-values for 3 methods at once: worst-rank wsr, worst-rank wmw, survivors only wsr
simulate_exp <- function(n, q2, HR, Delta, time, rho_X, rho_t){
  lambda2 = -1 / time * log(q2)
  lambda1 = lambda2 * HR
  #q1 = exp(-lambda1 * t)
  
  corr_t = matrix(c(1, rho_t, rho_t, 1), nrow = 2)
  t = rmvexp(n, rate=c(lambda1, lambda2), corr=corr_t)
  t1 = t[,1]
  t2 = t[,2]
  
  corr_X = matrix(c(1, rho_X, rho_X, 1), nrow = 2)
  X = rmvnorm(n, mean = c(0, 2^0.5 * Delta), cov = corr_X)
  X1 = X[,1]
  X2 = X[,2]
  
  eta = min(X) - 1 - time - 2 *(max(X)-min(X)) # in order to be worst
  distmin = min(abs(X1-X2))/time/10 # less than minimum difference between X1 and X2 divided by time
  
  X1_tilde = X1 * (t1 >= time) + (eta + t1*distmin) * (t1 < time)
  X2_tilde = X2 * (t2 >= time) + (eta + t2*distmin) * (t2 < time)
  
  X1_survivorOnly = X1[t1 >= time & t2 >= time]
  X2_survivorOnly = X2[t1 >= time & t2 >= time]
  
  if (length(X1_survivorOnly) == 0) {
    return (c(-1,-1,-1))
  }
  return (c(wilcox.test(X1_tilde, X2_tilde, alternative="less", paired=T)$p.value,
            wilcox.test(X1_tilde, X2_tilde, alternative="less", paired=F)$p.value,
            wilcox.test(X1_survivorOnly, X2_survivorOnly, alternative="less", paired=T)$p.value))
}

#simulate_probs_exp simulates probabilities for computing theoretical power(i.e., p1, p2, p3, p4 of Hetmmansperger(1984))
simulate_probs_exp <- function(n, q2, HR, Delta, time, rho_X, rho_t){
  lambda2 = -1 / time * log(q2)
  lambda1 = lambda2 * HR
  #q1 = exp(-lambda1 * t)
  
  corr_t = matrix(c(1, rho_t, rho_t, 1), nrow = 2)
  t = rmvexp(n, rate=c(lambda1, lambda2), corr=corr_t)
  t1 = t[,1]
  t2 = t[,2]
  
  corr_X = matrix(c(1, rho_X, rho_X, 1), nrow = 2)
  X = rmvnorm(n, mean = c(0, 2^0.5 * Delta), cov = corr_X)
  X1 = X[,1]
  X2 = X[,2]
  
  eta = min(X) - 1 - time - 2 *(max(X)-min(X)) # in order to be worst
  distmin = min(abs(X1-X2))/time/10 # less than minimum difference between X1 and X2 divided by time
  
  X1_tilde = X1 * (t1 >= time) + (eta + t1*distmin) * (t1 < time)
  X2_tilde = X2 * (t2 >= time) + (eta + t2*distmin) * (t2 < time)
  
  Y_tilde = X1_tilde - X2_tilde
  Y_tilde2 = permute(Y_tilde)
  Y_tilde3 = permute(Y_tilde)
  
  p1 = mean(Y_tilde > 0)
  p2 = mean((Y_tilde + Y_tilde2) > 0)
  p3 = mean(Y_tilde > 0 & ((Y_tilde + Y_tilde2) > 0)) # p3 = (p1 ** 2 + p2) / 2
  p4 = mean(((Y_tilde + Y_tilde2) > 0) * ((Y_tilde + Y_tilde3) > 0))
  return (c(p1,p2,p3,p4))
}



rmvnorm <- function(n, mean=NULL, cov=NULL) 
{
  ## munge parameters PRN and deal with the simplest univariate case
  
  if(is.null(mean))
    if(is.null(cov))
      return(rnorm(n))
  else
    mean = rep(0, nrow(cov))
  else if (is.null(cov))
    cov = diag(length(mean))
  
  ## gather statistics, do sanity checks
  
  D = length(mean)
  if (D != nrow(cov) || D != ncol(cov)) 
    stop("length of mean must equal nrow and ncol of cov")
  
  E = eigen(cov, symmetric=TRUE)
  if (any(E$val < 0)) 
    stop("Numerically negative definite covariance matrix")
  
  ## generate values and return
  
  mean.term = mean
  covariance.term = E$vec %*% (t(E$vec) * sqrt(E$val))
  independent.term = matrix(rnorm(n*D), nrow=D)
  
  drop(t(mean.term + (covariance.term %*% independent.term)))
}

## rmvexp: sampling from multivariate exponential distribution
rmvexp <- function(n, rate=1, corr=diag(length(rate)))
{
  ## extract parameters, do sanity checks, deal with univariate case
  
  if(!is.matrix(corr) || !isSymmetric(corr))
    stop("'corr' must be a symmetric matrix")
  D = ncol(corr)
  
  Dr = length(rate)
  if(Dr > D)
    warning("'rate' longer than width of 'corr', truncating to fit")
  if(Dr != D)
    rate = rep(rate, length.out=D)
  
  if(D == 1) rexp(n, rate)
  
  ## generate standard multivariate normal matrix, convert to CDF
  
  Z = rmvnorm(n, cov=corr)
  cdf = pnorm(Z)
  
  ## convert to exp, return
  
  sapply(1:D, function(d) qexp(cdf[,d], rate[d]))
}

rmvweibull <- function(n, shape=1, scale=1, corr=diag(length(shape)))
{
  ## extract parameters, do sanity checks, deal with univariate case
  
  if(!is.matrix(corr) || !isSymmetric(corr))
    stop("'corr' must be a symmetric matrix")
  D = ncol(corr)
  
  Dk = length(shape)
  if(Dk > D)
    warning("'shape' longer than width of 'corr', truncating to fit")
  if(Dk != D)
    shape = rep(shape, length.out=D)
  
  Db = length(scale)
  if(Db > D)
    warning("'scale' longer than width of 'corr', truncating to fit")
  if(Db != D)
    scale = rep(scale, length.out=D)
  
  if(D == 1) rweibull(n, shape, scale)
  
  ## generate standard multivariate normal matrix, convert to CDF

  Z = rmvnorm(n, cov=corr)
  cdf = pnorm(Z)
  
  ## convert to Weibull (WeiSD), return
  
  sapply(1:D, function(d) qweibull(cdf[,d], shape[d], scale[d]))
}

rmvloglogis <- function(n, shape=1, scale=1, corr=diag(length(shape))) 
{
  if(!is.matrix(corr) || !isSymmetric(corr))
    stop("'corr' must be a symmetric matrix")
  D = ncol(corr)
  
  Da = length(shape)
  if(Da > D)
    warning("'shape' longer than width of 'corr', truncating to fit")
  if(Da != D)
    shape = rep(shape, length.out=D)
  
  Db = length(scale)
  if(Db > D)
    warning("'scale' longer than width of 'corr', truncating to fit")
  if(Db != D)
    scale = rep(scale, length.out=D)
  
  if(D == 1) flexsurv::rllogis(n, shape, scale)
  
  ## generate standard multivariate normal matrix, convert to CDF
  
  Z = rmvnorm(n, cov=corr)
  cdf = pnorm(Z)
  
  ## convert to Weibull (WeiSD), return
  
  sapply(1:D, function(d) flexsurv::qllogis(cdf[,d], shape[d], scale[d]))
}

