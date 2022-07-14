library(sn)
library(nloptr)

### random sample from skew-normal copula ###
#############################################

rSNcopula=function(n,alpha,rho){
  
  # we assume Omega_bar=Omega
  
  # n: sample size
  
  ### multivariate sn parameters:
  # alpha: d vector with skewing parameters of multivariate skew-normal 
  # rho: d*(d-1)/2 vector with correlations of multivariate skew-normal
  #      ordered by column (rho_{2,1},rho_{3,1,..,rho_{d,1},rho_{3,2},..,rho_{d,d-1}
  
  # dimensions
  d=length(alpha)

  # construct correlation matrix of multivariate skew-normal
  Omega=diag(d)
  Omega[upper.tri(Omega)] <- Omega[lower.tri(Omega)] <- rho

  # # calculate delta from alpha and Omega
  # delta = Omega%*%alpha%*%(1+t(alpha)%*%Omega%*%alpha)^{-1/2}
  
  # generate sample from multivariate skew-normal
  x=rmsn(n = n,Omega = Omega,alpha = alpha)
  
  # calculation of marginal skewness parameters according to 
  # Azzalini & Capitanio 2014 chapter 5
  lambda=rep(NA,d)
  for(i in 1:d){
    # Omega1=1 because of diagonal element of correlation matrix
    Omega2=Omega[-i,-i]
    Omega12=matrix(Omega[i,-i],nrow=1)
    Omega221=Omega2-t(Omega12)%*%Omega12
    
    alpha1=alpha[i]
    alpha2=matrix(alpha[-i],ncol=1)
    
    lambda[i]=(1+t(alpha2)%*%Omega221%*%alpha2)^{-1/2}*(alpha1+Omega12%*%alpha2)
  }
  
  # calculate pseudo-observations
  u=matrix(NA,nrow=n,ncol=d)
  for(i in 1:d){
    u[,i]=psn(x = x[,i],omega = Omega[i,i],alpha = lambda[i])
  }
  
  return(u)
}




### Skew-normal copula distribution ###
#######################################

pSNcopula=function(u,alpha,rho,log=FALSE){
  # u: nxd matrix containing pseudo observations (in [0,1])
  # alpha: d vector containing multivariate skewing parameters
  # rho: d*(d-1) vector containing correlation of multivariate skew normal 
  #      distribution by column. 
  # log: should log density be returned or not
  
  n=nrow(u)
  d=ncol(u)
  
  # construct correlation matrix
  Omega=diag(d)
  Omega[upper.tri(Omega)] <- Omega[lower.tri(Omega)] <- rho
  
  # # calculate delta from alpha and Omega
  # delta = Omega%*%alpha%*%(1+t(alpha)%*%Omega%*%alpha)^{-1/2}
  
  # calculation of marginal skewness parameters according to 
  # Azzalini & Capitanio 2014 chapter 5
  lambda=rep(NA,d)
  for(i in 1:d){
    # Omega1=1 because of diagonal element of correlation matrix
    Omega2=Omega[-i,-i]
    Omega12=matrix(Omega[i,-i],nrow=1)
    Omega221=Omega2-t(Omega12)%*%Omega12
    
    alpha1=alpha[i]
    alpha2=matrix(alpha[-i],ncol=1)
    
    lambda[i]=(1+t(alpha2)%*%Omega221%*%alpha2)^{-1/2}*(alpha1+Omega12%*%alpha2)
  }
  
  # skew-normal quantile function evaluated in u (F^-1(u)) and
  # skew-normal distribution evaluated in skew-normal quantile in u (f(F^-1(u)))
  Finv=matrix(NA,nrow = n,ncol = d)
  for(i in 1:d){
    Finv[,i]=qsn(p = u[,i],omega = Omega[i,i],alpha = lambda[i])
  }
  
  # calculate density of multivariate skew-normal copula
  prob=pmsn(x = Finv,Omega = Omega,alpha = alpha)
  
  return(prob)
  
}





### Skew-normal copula density ###
##################################

dSNcopula=function(u,alpha,rho,log=FALSE){
  # u: nxd matrix containing pseudo observations (in [0,1])
  # alpha: d vector containing multivariate skewing parameters
  # rho: d*(d-1) vector containing correlation of multivariate skew normal 
  #      distribution by column. 
  # log: should log density be returned or not
  
  n=nrow(u)
  d=ncol(u)
  
  # construct correlation matrix
  Omega <- diag(x = 1/2,nrow = d)
  Omega[lower.tri(Omega)] <- rho
  Omega <- Omega + t(Omega)
  
  # # calculate delta from alpha and Omega
  # delta = Omega%*%alpha%*%(1+t(alpha)%*%Omega%*%alpha)^{-1/2}
  
  # calculation of marginal skewness parameters according to 
  # Azzalini & Capitanio 2014 chapter 5
  lambda=rep(NA,d)
  for(i in 1:d){
    # Omega1=1 because of diagonal element of correlation matrix
    Omega2=Omega[-i,-i]
    Omega12=matrix(Omega[i,-i],nrow=1)
    Omega221=Omega2-t(Omega12)%*%Omega12
    
    alpha1=alpha[i]
    alpha2=matrix(alpha[-i],ncol=1)
    
    lambda[i]=(1+t(alpha2)%*%Omega221%*%alpha2)^{-1/2}*(alpha1+Omega12%*%alpha2)
  }
  
  # skew-normal quantile function evaluated in u (F^-1(u)) and
  # skew-normal density evaluated in skew-normal quantile in u (f(F^-1(u)))
  Finv=matrix(NA,nrow = n,ncol = d)
  dn=matrix(NA,nrow = n,ncol = d)
  for(i in 1:d){
    Finv[,i]=qsn(p = u[,i],omega = Omega[i,i],alpha = lambda[i],engine="biv.nt.prob")
    dn[,i]=dsn(x = Finv[,i],omega = Omega[i,i],alpha = lambda[i])
  }
  
  # calculate density of multivariate skew-normal copula
  dens=dmsn(x = Finv,Omega = Omega,alpha = alpha)/apply(dn,1,prod)
  
  if(log==TRUE){
    return(log(dens))
  } else {
    return(dens)
  }
  
}




### Maximum-likelihood fit of skew-normal copula on pseudo-observations ###
###########################################################################

# negative log likelihood function
mll=function(param,u){
  d=ncol(u)
  alpha=param[1:d]
  rho=param[(d+1):(d*(d+1)/2)]
  logl=dSNcopula(u = u,alpha = alpha,rho = rho,log = TRUE)
  
  return(-sum(logl))
}

# negative log likelihood function for symmetric model
mllSym=function(rho,u){
  d=ncol(u)
  logl=dSNcopula(u = u,alpha = rep(0,d),rho = rho,log = TRUE)
  
  return(-sum(logl))
}

fitSNcopula=function(u,symmetric=FALSE,start=NULL,nstart=1,random=FALSE,...){

  # u: a nxd matrix containing the pseudo-observations (in [0,1])

  # symmetric: logical, should symmetry be assumed (alpha=0)

  # start: kx(d*(d+1)/2) matrix containing starting values for the parameters

  # nstart: if start is not provided, nstart random starting values are
  # generated in a plausible range. Defaults to 1

  # ...: options for the optimizer

  n=nrow(u)
  d=ncol(u)

  if(symmetric==TRUE){

    # starting values for minimizing negative log-likelihood
    if(is.null(start)){
      start=matrix(runif(n = nstart*(d*(d-1)/2),min = -0.25,max = 0.25),nrow = nstart,ncol=d*(d-1)/2)
    } else {
      nstart=nrow(start)
    }

    # lower and upper bounds for parameters
    lower=c(rep(-0.99,d*(d-1)/2))
    upper=c(rep(0.99,d*(d-1)/2))

    # minimizing the log-likelihood
    fits=vector(mode = "list", length = nstart)
    logl=rep(NA,nstart)
    for(i in 1:nstart){
      try({
        fit=optim(par = start[i,],fn = mllSym,method="L-BFGS-B",lower = lower,upper = upper,
                  u=u,...)
        fits[[i]]=fit
        logl[i]=fit$value
      },silent=T)
    }
  } else {
    # try to get good initial values not at random
    
    check=F
    try({
      inifit <- optim(par = runif(d*(d-1)/2,-1,1),fn = mllSym,method="L-BFGS-B",
                      lower = rep(-0.8,d*(d-1)/2),upper = rep(0.8,d*(d-1)/2),u=u,...)
      check=T
    },silent=T)
    if(is.null(start)){
      if(check==T & random==F){
        nstart=1
        alphastart=matrix(runif(n = nstart*d,min = -2,max = 2),nrow = nstart,ncol=d)
        rhostart=matrix(inifit$par,nrow=1,ncol=d*(d-1)/2)
      } else {
        alphastart=matrix(runif(n = nstart*d,min = -2,max = 2),nrow = nstart,ncol=d)
        rhostart=matrix(runif(n = nstart*d*(d-1)/2,min = -0.25,max = 0.25),nrow = nstart,ncol=d*(d-1)/2)
      }
      start=cbind(alphastart,rhostart)
    } else {
      nstart=nrow(start)
    }

    # lower and upper bounds for parameters
    lower=c(rep(-200,d),rep(-0.99,d*(d-1)/2))
    upper=c(rep(200,d),rep(0.99,d*(d-1)/2))

    # minimizing the log-likelihood
    fits=vector(mode = "list", length = nstart)
    logl=rep(NA,nstart)
    for(i in 1:nstart){
      try({
        fit=optim(par = start[i,],fn = mll,method="L-BFGS-B",lower = lower,upper = upper,
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