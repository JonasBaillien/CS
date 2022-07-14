library(ks)
fitNP=function(X,symmetric=TRUE,ub=NULL,lb=NULL,IF=2,GP=100){
  
  # fit a KDE to X over the range lb-ub
  # only works up to 6 dimensions
  
  # X: nxd matrix containing the observation in the rows
  # symmetric: logical, should symmetry be assumed or not
  # lb, ub: d-vector with lower bounds for the domain of the kde
  # IF: positive real number for the inflation of the bandwidth matrix of the KDE
  # GP: positive number for the gridpoints on one dimension for the KDE 
  
  # dimensions and data domain
  n <- nrow(X)
  d <- ncol(X)
  ranges <- apply(X,2,range)
  
  # checks and domain
  if(d>6){
    stop("Dimensions too high")
  }
  # 30% wider than original range of data 
  if(is.null(lb)){
    lb <- 1.3*ranges[1,]-0.3*ranges[2,] 
  }
  if(is.null(ub)){
    ub <- 1.3*ranges[2,]-0.3*ranges[1,]
  }
  
  # fitting the KDE
  if(symmetric==FALSE){
    fitReg <- kde(x = X,gridsize = GP,H = IF*Hpi(X),xmin = lb,xmax = ub,eval.points = X)
    logl <- sum(log(fitReg$estimate))
    return(list("dens"=fitReg$estimate,"bw"=fitReg$H,"ll"=logl))
  } else {
    fitPrim <- kde(x = X,gridsize = GP,H = IF*Hpi(X),xmin = lb,xmax = ub)
    theta.ind <- which(fitPrim$estimate == max(fitPrim$estimate),arr.ind = T)
    theta <- rep(NA,d)
    for(j in 1:d){
      theta[j] <- fitPrim$eval.points[[j]][theta.ind[j]]
    }

    Xsym <- sweep(x = -X,MARGIN = 2,STATS = 2*theta,FUN = "+")
    Xaug <- rbind(X,Xsym)
    
    fitSym <- kde(x = Xaug,H = IF*Hpi(X),xmin = lb,xmax = ub,eval.points = X)
    logl <- sum(log(fitSym$estimate))
    return(list("dens"=fitSym$estimate,"bw"=fitSym$H,"ll"=logl))
  }
}