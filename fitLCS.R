### function for fitting the linear combination of symmetric random variables ###
#################################################################################
source("~/PhD/code/lineaire combinatie/Code part 1 Functions.R")
fitLCS=function(data,basefunc,symmetric=FALSE,seed=NULL,maxiter=10^6,tol=10^-5,numstarts=20,start=NULL){
  # function for parameter estimation in a Multivariate Quantile Based Asymmetric Family of Distributions
  # model
  
  # data is a nxd-matrix which contains n-observations of d-variate points which form the data
  
  # basefunc is a vector  of length d which containing the names of univariate functions which 
  # are linearly combined to form the data. choises are "normal", "laplace", "logistic" or "t"
  
  # seed is a numeric seed to ensure consistent results on the same data
  
  # maxiter is a numeric value for the underlying optimizer which sets the maximum number of iterations
  # for the optimizer
  
  # tol is a numeric value for the underlying optimizer which sets the relative accuracy of the 
  # optimization step
  
  # numstarts is a numeric value which determines the number of different starting points for the optimizer
  
  # symmetric is a logical whether symmetry is assumed or not (this sets alll skewing parameters equal to 0.5)
  
  # start is an optional matrix of starting values for the parameter
  
  
  X=as.matrix(na.omit(data))
  
  # dimensions of X
  n=length(X[,1])
  d=length(X[1,])
  
  # check if correct number of base functions is supplied
  if(length(basefunc)!=d){
    stop("incorrect number of base functions supplied")
  }
  
  # check if basefunctions are ok
  possibilities=c("t","normal","logistic","laplace")
  test=!is.element(basefunc,possibilities)
  if(sum(test)>0){
    stop("incorrect basis function detected")
  }
  
  # generate starting values
  if(!is.null(seed)){
    set.seed(seed)
  } else {
    seed=sample(1:10^6,1)
    set.seed(seed)
  }
  
  # lower bounds for parameters in optimization
  lowermu=apply(X,2,min)
  lowerA=matrix(-sqrt(max(var(X))),nrow=d,ncol=d)
  lowerbounds=as.vector(cbind(lowermu,lowerA))
  
  # upper bounds for parameters in optimization
  uppermu=apply(X,2,max)
  upperA=matrix(sqrt(max(var(X))),nrow=d,ncol=d)
  upperbounds=as.vector(cbind(uppermu,upperA))
  
  if(is.null(start)){
    # for mu
    startmu=matrix(runif(d*numstarts,apply(X,2,min),apply(X,2,max)),nrow=numstarts,ncol=d,byrow=T)
    # for A 
    startA=matrix(runif(numstarts*d^2,min=-sqrt(max(var(X))),max=sqrt(max(var(X)))),nrow=numstarts,ncol=d^2,byrow=T)
  
    if(symmetric==T){
      alpha=rep(0.5,d)
    
      # combine starting values in matrix
      x0=cbind(startmu,startA)
      
    } else {
      loweralpha=rep(0.1,d)
      lowerbounds=as.vector(c(loweralpha,lowerbounds))
      upperalpha=rep(0.9,d)
      upperbounds=as.vector(c(upperalpha,upperbounds))
    
      # starting values for alpha
      startalpha=matrix(runif(d*numstarts,min = 0.2,max = 0.8),nrow=numstarts,ncol=d,byrow=T)
      
      # combine starting values in matrix
      x0=cbind(startalpha,startmu,startA)
    }
    
    if(is.element("t",basefunc)){
      
      # which basis functions are t-distributed
      indt=which(basefunc=="t")
      ldf=length(indt)
      
      # generate starting values for degrees of freedom as well as bounds
      lowerdf=rep(1,ldf)
      upperdf=rep(10000,ldf)
      startdf=matrix(runif(numstarts*ldf,min=2,max=100),ncol=ldf,nrow=numstarts,byrow=T)
      
      # pasting these to other starting values
      x0=cbind(x0,startdf)
      lowerbounds=c(lowerbounds,lowerdf)
      upperbounds=c(upperbounds,upperdf)
      
    } else {
      # in case no t-distribution is present
      indt=NULL
    }
  } else {
    
    x0=start

    if(symmetric==F){
      lowerbounds=as.vector(c(loweralpha,lowerbounds))
      upperbounds=as.vector(c(upperalpha,upperbounds))
    }
    
    if(is.element("t",basefunc)){
      
      # which basis functions are t-distributed
      indt=which(basefunc=="t")
      ldf=length(indt)
      
      # generate starting values for degrees of freedom as well as bounds
      lowerdf=rep(1,ldf)
      upperdf=rep(10000,ldf)
      
      lowerbounds=c(lowerbounds,lowerdf)
      upperbounds=c(upperbounds,upperdf)
      
    } else {
      # in case no t-distribution is present
      indt=NULL
    }
    
    checkL=sweep(x = start,MARGIN = 2,STATS = lowerbounds,FUN = "-")
    checkU=sweep(x = start,MARGIN = 2,STATS = upperbounds,FUN = "-")
    if(sum(checkL<0)>0 | sum(checkU>0)>0){
      stop("Starting values not in plausible range")
    }
  }
  
  M=ncol(x0)
  if(length(upperbounds)!=M){
    stop("Incorrect number of parameter values provided in starting values")
  }
  
  
  # holding matrix for best found parameters (one row is one set)
  paramsfit=data.frame(matrix(NA,ncol=M,nrow=numstarts))
  a=intToUtf8(945)
  m=intToUtf8(956) # character for mu
  if(symmetric==T){
    if(is.null(indt)){
      colnames(paramsfit)=c(paste0(m,1:d),paste0("A",apply(expand.grid(1:d, 1:d), 1, paste, collapse="")))
    } else {
      colnames(paramsfit)=c(paste0(m,1:d),paste0("A",apply(expand.grid(1:d, 1:d), 1, paste, collapse="")),paste0("df",indt))
    }
  } else {
    if(is.null(indt)){
      colnames(paramsfit)=c(paste0(a,1:d),paste0(m,1:d),paste0("A",apply(expand.grid(1:d, 1:d), 1, paste, collapse="")))
    } else {
      colnames(paramsfit)=c(paste0(a,1:d),paste0(m,1:d),paste0("A",apply(expand.grid(1:d, 1:d), 1, paste, collapse="")),paste0("df",indt))
    }
  }
  
  
  
  # vector with log-likelihood in found optimum 
  loglfit=rep(NA,numstarts)
  
  
  # main loop for parameter estimation with symmetry assumed
  for(i in 1:numstarts){
      
    try({
      # set seed
      set.seed(seed+i)
      # starting values for parameters
      parstart=x0[i,]
        
      # optimization
      output=bobyqa(x0 = parstart,fn = opt_funLCS,lower=lowerbounds,upper=upperbounds,X=X,basefunc=basefunc,indt=indt,symmetric=symmetric,control = list(maxeval=120000))
      paramsfit[i,]=output$par
      loglfit[i]=output$value
    },silent=T)
  }
  
  # only return starting value which provides best fit (lowest minus log-likelihood)
  indmin=which.min(loglfit)
  loglfit=-loglfit[indmin]
  
  if(symmetric==T){  
    paramsfit=c(rep(0.5,d),unlist(paramsfit[indmin,]))
    names(paramsfit)[1:d]=paste0(a,1:d)
  } else {
    paramsfit=paramsfit[indmin,]
  }
  
  return(list("fitted parameters"=paramsfit,"log likelihood fit"=loglfit))
}



# minus log-likelihood function which is to be minimized for symmetric components
opt_funLCS=function(parms,X,basefunc,symmetric=TRUE,indt=NULL){
  # parms is a vector containing the parameter values. In order these are: alpha, mu, A (by column) and 
  # possibly degrees of freedom for student-t distributions. 
  
  # X is a nxd matrix containing the data
  
  # basefunc is a vector  of length d which containing the names of univariate functions which 
  # are linearly combined to form the data. choices are "normal", "laplace", "logistic" or "t". 
  # indices of t-distributions within this vector are supplied via indt. This is used to generate
  # a vector of length d which contains the degrees of freedom with the corresponding t-distribution in
  # in basefunc
  
  # if symmetric=TRUE, alpha is not estimated and assumed to be 0.5
  # otherwise, alpha is contained as the first d elements of parms
  
  # dimensions
  n=length(X[,1])
  d=length(X[1,])
  
  # parameter values
  if(symmetric==TRUE){
    # regular parameters
    alpha=rep(0.5,d)
    mu=parms[1:d]
    A=matrix(parms[(d+1):(d^2+d)],nrow=d)
    
    # degrees of freedom
    if(!is.null(indt)){
      tpars=rep(NA,d)
      tpars[indt]=parms[(d^2+d+1):(d^2+d+length(indt))]
    }
  } else {
    # regular parameters
    alpha=parms[1:d]
    mu=parms[(d+1):(2*d)]
    A=matrix(parms[(2*d+1):(d^2+2*d)],nrow=d)
    
    # creating vector with degrees of freedom
    if(!is.null(indt)){
      tpars=rep(NA,d)
      tpars[indt]=parms[(d^2+2*d+1):(d^2+2*d+length(indt))]
    }
  }
  B=solve(A)
  
  
  
  # X*A^{-1}
  V=X%*%B
  # mu*A^{-1}
  VV=mu%*%B
  
  # calculation of minus log-likelihood, max ensures that no inf's of NA's are produced
  value=-n*log(abs(det(B)))
  for(i in 1:d){
    if(basefunc[i]=="laplace"){
      value=value-sum(log(twopiecelaplacedensity(V[,i]-VV[i],alpha[i])))
    } else if(basefunc[i]=="normal"){
      value=value-sum(log(twopiecenormaldensity(V[,i]-VV[i],alpha[i])))
    } else if(basefunc[i]=="logistic"){
      value=value-sum(log(twopiecelogisticdensity(V[,i]-VV[i],alpha[i])))
    } else {
      value=value-sum(log(twopiecestudentdensity(V[,i]-VV[i],alpha[i],nu=tpars[i])))
    }
  }
  
  return(value)
}