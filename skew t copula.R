library(sn)
library(nloptr)

### random sample from skew-t copula ###
########################################

rSTcopula=function(n,alpha,rho,nu){
  
  # we assume Omega_bar=Omega
  
  # n: sample size
  
  ### multivariate sn parameters:
  # alpha: d vector with skewing parameters of multivariate skew-normal 
  # rho: d*(d-1)/2 vector with correlations of multivariate skew-normal
  #      ordered by column (rho_{2,1},rho_{3,1,..,rho_{d,1},rho_{3,2},..,rho_{d,d-1}
  # nu: degrees of freedom of t-distribution; >d, numeric
  
  # dimensions
  d=length(alpha)
  
  # construct correlation matrix
  Omega = diag(x = 1/2,nrow = d)
  Omega[lower.tri(Omega)] = rho
  Omega = Omega + t(Omega)
  
  # calculate delta from alpha and Omega
  delta = Omega%*%alpha%*%(1+t(alpha)%*%Omega%*%alpha)^{-1/2}
  
  # generate sample from multivariate skew-t
  x=rmst(n = n,Omega = Omega,alpha = alpha,nu = nu)
  
  # calculation of marginal skewness parameters according to 
  # Azzalini & Capitanio 2014 chapter 6 & Yoshiba 2018
  zeta=delta/(sqrt(1-delta^2))
  
  # calculate pseudo-observations
  u=matrix(NA,nrow=n,ncol=d)
  for(i in 1:d){
    u[,i]=pst(x = x[,i],omega = Omega[i,i],alpha = zeta[i],nu = nu)
  }
  
  return(u)
}




### Skew-t copula distribution ###
##################################

pSTcopula=function(u,alpha,rho,nu){
  # u: nxd matrix containing pseudo observations (in [0,1])
  # alpha: d vector containing multivariate skewing parameters
  # rho: d*(d-1) vector containing correlation of multivariate skew normal 
  #      distribution by column. 
  # nu: degrees of freedom of t-distribution >d numeric
  
  n=nrow(u)
  d=ncol(u)
  
  # construct correlation matrix 
  Omega = diag(x = 1/2,nrow = d)
  Omega[lower.tri(Omega)] = rho
  Omega = Omega + t(Omega)
  
  # calculate delta from alpha and Omega
  delta = Omega%*%alpha%*%(1+t(alpha)%*%Omega%*%alpha)^{-1/2}
  
  # calculation of marginal skewness parameters according to 
  # Azzalini & Capitanio 2014 chapter 6 & Yoshiba 2018
  zeta=delta/(sqrt(1-delta^2))
  
  
  # skew-normal quantile function evaluated in u (F^-1(u)) and
  # skew-normal distribution evaluated in skew-normal quantile in u (f(F^-1(u)))
  Finv=matrix(NA,nrow = n,ncol = d)
  for(i in 1:d){
    Finv[,i]=qst(p = u[,i],omega =1,alpha = zeta[i],nu = nu)
  }
  
  # calculate density of multivariate skew-t copula 
  # pmst only accepts a vector, not a matrix
  prob=apply(Finv,1,pmst,Omega = Omega,alpha = alpha,nu = nu)
  
  return(prob)
  
}





### Skew-normal copula density ###
##################################

dSTcopula=function(u,alpha,rho,nu,log=FALSE){
  # u: nxd matrix containing pseudo observations (in [0,1])
  # alpha: d vector containing multivariate skewing parameters
  # rho: d*(d-1) vector containing correlation of multivariate skew normal 
  #      distribution by column. 
  # nu: degrees of freedom of t-distribution >d numeric
  # log: should log density be returned or not
  
  n=nrow(u)
  d=ncol(u)
  
  # construct correlation matrix 
  Omega = diag(x = 1/2,nrow = d)
  Omega[lower.tri(Omega)] = rho
  Omega = Omega + t(Omega)
  
  # calculate delta from alpha and Omega
  delta = Omega%*%alpha%*%(1+t(alpha)%*%Omega%*%alpha)^{-1/2}
  
  # calculation of marginal skewness parameters according to 
  # Azzalini & Capitanio 2014 chapter 6 & Yoshiba 2018
  zeta=delta/(sqrt(1-delta^2))
  
  # skew-t quantile function evaluated in u (F^-1(u)) and
  # skew-t density evaluated in skew-t quantile in u (f(F^-1(u)))
  Finv=matrix(NA,nrow = n,ncol = d)
  dt=matrix(NA,nrow = n,ncol = d)
  for(i in 1:d){
    Finv[,i]=qst(p = u[,i],omega = 1,alpha = zeta[i],nu = nu)
    dt[,i]=dst(x = Finv[,i],omega = 1,alpha = zeta[i],nu = nu)
  }
  
  # calculate density of multivariate skew-normal copula
  dens=dmst(x = Finv,Omega = Omega,alpha = alpha,nu = nu)/apply(dt,1,prod)
  
  if(log==TRUE){
    return(log(dens))
  } else {
    return(dens)
  }
  
}




### Maximum-likelihood fit of skew-normal copula on pseudo-observations ###
###########################################################################

# negative log likelihood function
mllT=function(param,u){
  d=ncol(u)
  alpha=param[1:d]
  rho=param[(d+1):(d*(d+1)/2)]
  nu=param[(d*(d+1)/2)+1]
  logl=dSTcopula(u = u,alpha = alpha,rho = rho,nu = nu,log = TRUE)
  
  return(-sum(logl))
}

# negative log likelihood function for symmetric model
mllTSym=function(param,u){
  d=ncol(u)
  rho=param[1:(d*(d-1)/2)]
  nu=param[d*(d-1)/2+1]
  logl=dSTcopula(u = u,alpha = rep(0,d),rho = rho,nu = nu,log = TRUE)
  
  return(-sum(logl))
}

fitSTcopula=function(u,symmetric=FALSE,start=NULL,nstart=1,...){
  
  # u: a nxd matrix containing the pseudo-observations (in [0,1])
  
  # symmetric: logical, should symmetry be assumed (alpha=0)
  
  # start: kx(d*(d+1)/2+1) matrix containing starting values for the parameters;
  #        kx(d*(d-1)/2+1) matrix when symmetric=FALSE
  
  # nstart: if start is not provided, nstart random starting values are
  # generated in a plausible range. Defaults to 1
  
  # ...: options for the optimizer
  
  n=nrow(u)
  d=ncol(u)
  
  if(symmetric==TRUE){
    
    # starting values for minimizing negative log-likelihood
    if(is.null(start)){
      startrho=matrix(runif(n = nstart*(d*(d-1)/2),min = -0.55,max = 0.55),nrow = nstart,ncol=d*(d-1)/2)
      startnu=matrix(runif(n = nstart,min = d,max = 40),nrow=nstart,ncol=1)
      start=cbind(startrho,startnu)
    } else {
      nstart=nrow(start)
    }
    
    # lower and upper bounds for parameters
    lower=c(rep(-0.8,d*(d-1)/2),d)
    upper=c(rep(0.8,d*(d-1)/2),1000)
    
    # minimizing the log-likelihood
    fits=vector(mode = "list", length = nstart)
    logl=rep(NA,nstart)
    for(i in 1:nstart){
      try({
        fit=optim(par = start[i,],fn = mllTSym,method="L-BFGS-B",lower = lower,upper = upper,
                  u=u,...)
        fits[[i]]=fit
        logl[i]=fit$value
      },silent=T)
    }
  } else {
    # starting values for minimizing negative log-likelihood
    if(is.null(start)){
      alphastart=matrix(runif(n = nstart*d,min = -5,max = 5),nrow = nstart,ncol=d)
      rhostart=matrix(runif(n = nstart*(d*(d-1)/2),min = -0.55,max = 0.55),nrow = nstart,ncol=d*(d-1)/2)
      nustart=matrix(runif(n = nstart,min = d,max = 40),nrow=nstart,ncol=1)
      start=cbind(alphastart,rhostart,nustart)
    } else {
      nstart=nrow(start)
    }
    
    # lower and upper bounds for parameters
    lower=c(rep(-20,d),rep(-0.6,d*(d-1)/2),d)
    upper=c(rep(20,d),rep(0.6,d*(d-1)/2),60)
    
    # minimizing the log-likelihood
    fits=vector(mode = "list", length = nstart)
    logl=rep(NA,nstart)
    for(i in 1:nstart){
      try({
        fit=optim(par = start[i,],fn = mllT,method="L-BFGS-B",lower = lower,upper = upper,
                  u=u,...)
        fits[[i]]=fit
        logl[i]=fit$value
      },silent=T)
    }
  }
  
  index=which.min(logl)
  optimum=fits[[index]]$par
  return(list("par"=optimum,"ll"=-logl[index]))
  
}

# n=100
# alpha = c(5,0,5,0)
# rho = c(0.2,0.35,-0.1,0.15,0.2,-0.4)
# nu = 6
# X=rSTcopula(n = n,alpha = alpha,rho = rho,nu = nu)
# 
# dSTcopula(u = X,alpha = alpha,rho = rho,nu = nu,log = F)
# pSTcopula(u = X,alpha = alpha,rho = rho,nu = nu)
# fitSTcopula(u = X,symmetric = T,nstart = 2)
# fitSTcopula(u = X,symmetric = F,nstart = 2)
