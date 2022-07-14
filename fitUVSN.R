# function for fitting univariate skew normal distribution to data used in apply
fitUVSN=function(Y,symmetric=T){
  n <- length(Y)
  if(symmetric==T){
    fit <- selm.fit(x = matrix(1,nrow=n,ncol=1),y = Y,family = 'SN',fixed.param = list(alpha=0),selm.control = list(method="MLE"))
    return(list("par"=fit$param$dp.complete,"ll"=fit$logL))
  } else {
    fit <- selm.fit(x = matrix(1,nrow=n,ncol=1),y = Y,family = 'SN',selm.control = list(method="MLE"))
    return(list("par"=fit$param$dp,"ll"=fit$logL))
  }
}