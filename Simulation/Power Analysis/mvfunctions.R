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

