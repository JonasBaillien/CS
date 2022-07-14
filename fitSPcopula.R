### file  for determining the GLRT for symmetry in          ###   
### symmetric copula with non-parametric margins            ###
###############################################################
library(copula)
fitSPcopula=function(X,cop,theta=NULL,theta.method="median",bw=NULL,symmetric=FALSE){
  
  # function for semi-parametric fit of a (symmetric) copula to data
  
  # X is a nxd-matrix which contains n-observations of d-variate points which form the data
  
  # cop is an object of the copula-class. Since possible symmetry is required, options are
  # Frank (d=2) and Farlie-Gumbel-Morgenstern (d=2) elliptical copulas
  
  # theta is a d-vector containing the points of symmetry. if NULL, theta.method 
  # is used to find an estimate. Default is "median", other options are "mode" 
  # (in case unimodal) or "mean" (less robust)
  
  # bw is a d-vector containing the bandwidth for each margin. Only used when 
  # theta.method="mode" Default is NULL, in that case the nrd0-bandwidth is used
  
  # symmetric is a logical whether symmetry is assumed or not (this only has effect
  # on the estimation of the margins) if TRUE, the pseudo observations are generated
  # as the CDF they have if the data and the reflected data are combined (and thus the 
  # data itself is symmetric around theta)

  X=as.matrix(na.omit(X))
  
  n=nrow(X)
  d=ncol(X)
  
  # check for correct copula object
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
  
  # determine point of symmetry theta
  if(is.null(theta)){
    if(is.element(theta.method,c("median","mean","mode"))){
      if(theta.method=="median"){
        theta=apply(X,2,median)
      } else if(theta.method=="mean"){
        theta=apply(X,2,mean)
      } else {
  
        # set bandwidth vector
        if(is.null(bw)){
          bw=apply(X,2,nrd0)
        } else {
          if(length(bw)!=d){
            stop("Incorrect length of bandwidth matrix provided")
          }
        }
        
        theta=rep(NA,d)
        for(i in 1:d){
          y=density(x = X[,i],bw=bw[i])
          theta[i]=y$x[which.max(y$y)]
        }

      }
    } else {
      stop("Invalid theta.method provided")
    }
    
  } else {
    if(length(theta)!=d){
      stop("Length of theta does not match dimension")
    }
  }
  
  # generating pseudo observations
  if(symmetric==FALSE){
    U=pobs(X)
  } else {
    Xm=sweep(x = -X,MARGIN = 2,STATS = 2*theta,FUN = "+")
    U=pobs(rbind(X,Xm))[1:n,]
  }
  
  # fitting copula to pseudo observations
  fit=fitCopula(copula = cop,data = U)
  
  
  return(fit)
}