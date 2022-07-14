### file  for determining the GLRT for symmetry in          ###   
### symmetric copula with non-parametric margins            ###
###############################################################
library(copula)
source("C:/Users/u0125240/Documents/PhD/code/copulas/univariatefit.R")


fitCopulaQBAM=function(X,cop,margfunc=NULL,symmetric=FALSE){
  
  # function for parametric fit of a (symmetric) copula to data where margins 
  # are assumed to be from a QBA distribution
  
  # X is a nxd-matrix which contains n-observations of d-variate points which form the data
  
  # cop is an object of the copula-class. Since possible symmetry is required, options are
  # Frank (d=2) and Farlie-Gumbel-Morgenstern (d=2) elliptical copulas
  
  # margfunc is a d-vector containing the name of the QBA margin. Options are 
  # "normal", "laplace", "logistic" and "t". If NULL, the best fitting one based
  # on AIC is used from the four possibilities
  
  # symmetric is a logical whether symmetry is assumed or not (this only has effect
  # on the estimation of the margins) if TRUE, in the estimation of the QBA-margins
  # the skewness parameter is fixed at 0.5 (indicating symmetry) otherwise it is  
  # estimated based on the data provided
  
  X=as.matrix(na.omit(X))
  
  n=nrow(X)
  d=ncol(X)
  
  ### check for correct copula object
  check.cop=
    if(d==2){
      if(class(cop)!="frankCopula" & !is(cop,"ellipCopula") & class(cop)!="fgmCopula"){
        stop("Invalid copula object provided")
      }
    } else {
      if(!is(cop,"ellipCopula")){
        stop("Invalid copula object provided")
      }
    }
  

  
  
  ### estimating margins and pseudo-observations
  if(is.null(margfunc)){
    
    # holding matrices
    U=matrix(NA,nrow=n,ncol=d)
    
    parameters=matrix(NA,nrow=4,ncol=d)
    rownames(parameters)=c('alpha','mu','phi','nu')
    
    margin=rep(NA,d)
    
    loglikelihood=rep(NA,d)
    
    # choices for margin functions
    mf=c("normal","laplace","logistic","t")
    
    for(i in 1:d){
      # fitting the the possible margins
      fitN=fitSAND(X = X[,i],symmetric=symmetric)
      fitLa=fitSALaD(X = X[,i],symmetric=symmetric)
      fitLo=fitSALoD(X = X[,i],symmetric=symmetric)
      fitT=fitSATD(X = X[,i],symmetric=symmetric)
      
      # uniform tranform (CDF values) of the data
      Ut=matrix(NA,nrow=n,ncol=4)
      Ut[,1]=QBAsyDist::pAND(q=X[,i],mu=fitN$mu,phi=fitN$phi,alpha=fitN$alpha)
      Ut[,2]=QBAsyDist::pALaD(q=X[,i],mu=fitLa$mu,phi=fitLa$phi,alpha=fitLa$alpha)
      Ut[,3]=QBAsyDist::pALoD(q=X[,i],mu=fitLo$mu,phi=fitLo$phi,alpha=fitLo$alpha)
      Ut[,4]=QBAsyDist::pATD(q=X[,i],mu=fitT$mu,phi=fitT$phi,alpha=fitT$alpha,nu = fitT$nu)
      
      # fitted parameters
      pars=matrix(NA,nrow=4,ncol=4)
      rownames(pars)=c('alpha','mu','phi','nu')
      pars[1:3,1]=c(fitN$alpha,fitN$mu,fitN$phi)
      pars[1:3,2]=c(fitLa$alpha,fitLa$mu,fitLa$phi)
      pars[1:3,3]=c(fitLo$alpha,fitLo$mu,fitLo$phi)
      pars[1:4,4]=c(fitT$alpha,fitT$mu,fitT$phi,fitT$nu)
      
      # log-likelihood
      logl=rep(NA,4)
      logl[1]=fitN$LogLikelihood
      logl[2]=fitLa$LogLikelihood
      logl[3]=fitLo$LogLikelihood
      logl[4]=fitT$LogLikelihood
      
      AICfit=-2*c(fitN$LogLikelihood,fitLa$LogLikelihood,fitLo$LogLikelihood,fitT$LogLikelihood)+2*c(3,3,3,4)
      margin[i]=mf[which.min(AICfit)]
      U[,i]=Ut[,which.min(AICfit)]
      parameters[,i]=pars[,which.min(AICfit)]
      loglikelihood[i]=logl[which.min(AICfit)]
          
    }
  } else {
    
    # holding matrices
    U=matrix(NA,nrow=n,ncol=d)
      
    parameters=matrix(NA,nrow=4,ncol=d)
    rownames(parameters)=c('alpha','mu','phi','nu')
    
    margin=margfunc
    
    loglikelihood=rep(NA,d)
    
    for(i in 1:d){
      fitt=switch(margfunc[i],"normal"=fitSAND(X = X[,i],symmetric=symmetric),
                  "logistic"=fitSALoD(X = X[,i],symmetric=symmetric),
                  "laplace"=fitSALaD(X = X[,i],symmetric=symmetric),
                  "t"=fitSATD(X = X[,i],symmetric=symmetric))
      U[,i]=switch(margfunc[i],"normal"=QBAsyDist::pAND(q=X[,i],mu=fitt$mu,phi=fitt$phi,alpha=fitt$alpha),
                        "logistic"=QBAsyDist::pALoD(q=X[,i],mu=fitt$mu,phi=fitt$phi,alpha=fitt$alpha),
                        "laplace"=QBAsyDist::pALaD(q=X[,i],mu=fitt$mu,phi=fitt$phi,alpha=fitt$alpha),
                        "t"=QBAsyDist::pATD(q=X[,i],mu=fitt$mu,phi=fitt$phi,alpha=fitt$alpha,nu = fitt$nu))
      parameters[1:3,i]=c(fitt$alpha,fitt$mu,fitt$phi)
      if(margfunc[i]=="t"){
        parameters[4,i]=fitt$nu
      }
      loglikelihood[i]=fitt$LogLikelihood
    }
    
  }
  
  fitM=list("margins"=margin,"parameters"=parameters,"LogLikelihood"=loglikelihood)
  
  # fitting copula to pseudo observations
  fitC=fitCopula(copula = cop,data = U)
  
  
  return(list("fitC"=fitC,"fitM"=fitM))
}


### normal fitting
fitSAND=function(X,symmetric=FALSE){
  # function for fitting a quantile-based normal distribution to the data
  # using maximum likelihood
  
  # X is a numeric vector containing the data
  # symmetric is a logical indicating if symmetry is assumed or not (alpha=0.5)
    
  # upper and lower bounds for parameters in optimization (alpha,mu,phi)
  if(symmetric==FALSE){
    lowerbounds=c(0,-Inf,0)
    upperbounds=c(1,Inf,Inf)
  } else {
    lowerbounds=c(0.5,-Inf,0)
    upperbounds=c(0.5,Inf,Inf)
  }
    
  # generate starting values
  if(symmetric==FALSE){
    startalpha=runif(10)
  } else {
    startalpha=rep(0.5,10)
  }
  startmu=runif(10,min=min(X),max=max(X))
  startphi=runif(10,min=0,max=sd(X))
  
  # combine starting values in matrix
  x0=cbind(startalpha,startmu,startphi)
  
  # holding vectors for parameter estimates and log-likelihood
  parsfit=matrix(NA,nrow=10,ncol=3)
  loglfit=rep(NA,10)
      
  # main loop for parameter estimation using the bobyqa function from the nloptr package
  for(i in 1:10){
        
      # in case of error for certain starting values, they are suppressed
      try({
        # minimization of minus the log likelihood
        output=nloptr::bobyqa(x0=x0[i,],fn=dQBN,lower=lowerbounds,upper=upperbounds,nl.info=F,x=X,control = list(maxeval=50000,xtol_rel=10^-5))
        i=i+1
        
        # optimal parameters and minus log-likelihood
        parsfit[i,]=output$par
        loglfit[i]=output$value
      },silent=T)
    }
      
  # returning best fit
  indmin=which.min(loglfit)
  bestpars=parsfit[indmin,]
      
  return(list("alpha"=bestpars[1],"mu"=bestpars[2],"phi"=bestpars[3],"LogLikelihood"=-loglfit[indmin]))
}


fitSALaD=function(X,symmetric=FALSE){
  # function for fitting a quantile-based Laplace distribution to the data
  # using maximum likelihood
  
  # X is a numeric vector containing the data
  # symmetric is a logical indicating if symmetry is assumed or not (alpha=0.5)
  
  
  # upper and lower bounds for parameters in optimization (alpha,mu,phi)
  if(symmetric==FALSE){
    lowerbounds=c(0,-Inf,0)
    upperbounds=c(1,Inf,Inf)
  } else {
    lowerbounds=c(0.5,-Inf,0)
    upperbounds=c(0.5,Inf,Inf)
  }
  
  # generate starting values
  if(symmetric==FALSE){
    startalpha=runif(10)
  } else {
    startalpha=rep(0.5,10)
  }
  startmu=runif(10,min=min(X),max=max(X))
  startphi=runif(10,min=0,max=sd(X))
  
  # combine starting values in matrix
  x0=cbind(startalpha,startmu,startphi)
  
  # holding vectors for parameter estimates and log-likelihood
  parsfit=matrix(NA,nrow=10,ncol=3)
  loglfit=rep(NA,10)
  
  # main loop for parameter estimation using the bobyqa function from the nloptr package
  for(i in 1:10){
      
    # in case of error for certain starting values, they are suppressed
    try({

      # minimization of minus the log likelihood
      output=bobyqa(x0=x0[i,],fn=dQBLa,lower=lowerbounds,upper=upperbounds,nl.info=F,x=X,control = list(maxeval=50000,xtol_rel=10^-5))
      
      # optimal parameters and minus log-likelihood
      parsfit[i,]=output$par
      loglfit[i]=output$value
    },silent=T)
  }
    
  # returning best fit
  indmin=which.min(loglfit)
  bestpars=parsfit[indmin,]
    
  return(list("alpha"=bestpars[1],"mu"=bestpars[2],"phi"=bestpars[3],"LogLikelihood"=-loglfit[indmin]))
}
  
  
fitSALoD=function(X,symmetric=FALSE){
  # function for fitting a quantile-based logistic distribution to the data
  # using maximum likelihood
  
  # X is a numeric vector containing the data
  # symmetric is a logical indicating if symmetry is assumed or not (alpha=0.5)
  
  
  # upper and lower bounds for parameters in optimization (alpha,mu,phi)
  if(symmetric==FALSE){
    lowerbounds=c(0,-Inf,0)
    upperbounds=c(1,Inf,Inf)
  } else {
    lowerbounds=c(0.5,-Inf,0)
    upperbounds=c(0.5,Inf,Inf)
  }
  
  # generate starting values
  if(symmetric==FALSE){
    startalpha=runif(10)
  } else {
    startalpha=rep(0.5,10)
  }
  startmu=runif(10,min=min(X),max=max(X))
  startphi=runif(10,min=0,max=sd(X))
  
  # combine starting values in matrix
  x0=cbind(startalpha,startmu,startphi)
  
  # holding vectors for parameter estimates and log-likelihood
  parsfit=matrix(NA,nrow=10,ncol=3)
  loglfit=rep(NA,10)
  
  # main loop for parameter estimation using the bobyqa function from the nloptr package
  for(i in 1:10){
    
    # in case of error for certain starting values, they are suppressed
    try({
      
      # minimization of minus the log likelihood
      output=bobyqa(x0=x0[i,],fn=dQBLo,lower=lowerbounds,upper=upperbounds,nl.info=F,x=X,control = list(maxeval=50000,xtol_rel=10^-5))
      
      # optimal parameters and minus log-likelihood
      parsfit[i,]=output$par
      loglfit[i]=output$value
    },silent=T)
  }
  
  # returning best fit
  indmin=which.min(loglfit)
  bestpars=parsfit[indmin,]
  
  return(list("alpha"=bestpars[1],"mu"=bestpars[2],"phi"=bestpars[3],"LogLikelihood"=-loglfit[indmin]))
}


fitSATD=function(X,symmetric=FALSE){
  # function for fitting a quantile-based student's t-distribution to the data
  # using maximum likelihood
  
  # X is a numeric vector containing the data
  # symmetric is a logical indicating if symmetry is assumed or not (alpha=0.5)
  
  
  # upper and lower bounds for parameters in optimization (alpha,mu,phi)
  if(symmetric==FALSE){
    lowerbounds=c(0,-Inf,0,0)
    upperbounds=c(1,Inf,Inf,1000)
  } else {
    lowerbounds=c(0.5,-Inf,0,0)
    upperbounds=c(0.5,Inf,Inf,1000)
  }
  
  
  
  # generate starting values
  if(symmetric==FALSE){
    startalpha=runif(10)
  } else {
    startalpha=rep(0.5,10)
  }
  startmu=runif(10,min=min(X),max=max(X))
  startphi=runif(10,min=0,max=sd(X))
  startnu=runif(10,min=2,max=200)
  
  # combine starting values in matrix
  x0=cbind(startalpha,startmu,startphi,startnu)
  
  # holding vectors for parameter estimates and log-likelihood
  parsfit=matrix(NA,nrow=10,ncol=4)
  loglfit=rep(NA,10)
  
  # main loop for parameter estimation using the bobyqa function from the nloptr package
  for(i in 1:10){
    
    # in case of error for certain starting values, they are suppressed
    try({
      
      # minimization of minus the log likelihood
      output=bobyqa(x0=x0[i,],fn=dQBT,lower=lowerbounds,upper=upperbounds,nl.info=F,x=X,control = list(maxeval=50000,xtol_rel=10^-5))
      
      # optimal parameters and minus log-likelihood
      parsfit[i,]=output$par
      loglfit[i]=output$value
    },silent=T)
  }
  
  # returning best fit
  indmin=which.min(loglfit)
  bestpars=parsfit[indmin,]
  
  return(list("alpha"=bestpars[1],"mu"=bestpars[2],"phi"=bestpars[3],"nu"=bestpars[4],"LogLikelihood"=-loglfit[indmin]))
}




