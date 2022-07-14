### function for determining the best fitting quantile- ###
### based margin to fit a margin of the data            ###
###########################################################

library(QBAsyDist)
source(file="~/PhD/code/copulas/univariatefit.R")
source(file="C:/Users/u0125240/Documents/PhD/code/symmetry/fitPcopula.R")

marginfit=function(data,crit="AIC",all=F,symmetric=F){
  # accepts a numeric vector as input which forms a single margin of the data on which 
  # the model is fit. The best suited margin is determined from a list of 4 and based
  # on either the correlation of a uniform QQ-plot (crit="QQ") of the fitted distribution
  # or AIC/BIC (crit="AIC"/"BIC") or the best one according to the p-value of a Kolmogorov-
  # Smirnov test. If all=T the criterion for all possibilities is returned so the user can 
  # choose on his/her preference, otherwise only the best one toghether with the uniform 
  # transformation of the data, the fitted parameters and the log-likelihood is returned
  
  
  # length of data vector
  n=length(data)
  
  # choices for margin functions
  mf=c("normal","laplace","logistic","t")
  
  # fitting the the possible margins
  fitN=fitSAND(data,symmetric=symmetric)
  fitLa=fitSALaD(data,symmetric=symmetric)
  fitLo=fitSALoD(data,symmetric=symmetric)
  fitT=fitSATD(data,symmetric=symmetric)
  
  # uniform tranform (CDF values) of the data
  U=matrix(NA,nrow=n,ncol=4)
  U[,1]=QBAsyDist::pAND(q=data,mu=fitN$mu,phi=fitN$phi,alpha=fitN$alpha)
  U[,2]=QBAsyDist::pALaD(q=data,mu=fitLa$mu,phi=fitLa$phi,alpha=fitLa$alpha)
  U[,3]=QBAsyDist::pALoD(q=data,mu=fitLo$mu,phi=fitLo$phi,alpha=fitLo$alpha)
  U[,4]=QBAsyDist::pATD(q=data,mu=fitT$mu,phi=fitT$phi,alpha=fitT$alpha,nu = fitT$nu)
  
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
  
  if(crit=="AIC"){
    
    AICfit=-2*c(fitN$LogLikelihood,fitLa$LogLikelihood,fitLo$LogLikelihood,fitT$LogLikelihood)+2*c(3,3,3,4)
    if(all==T){
      
      return(cbind("margin"=mf,"AIC"=AICfit,"parameters"=t(pars)))
      
    } else {
      
      return(list("margin"=mf[which.min(AICfit)],
                  "U"=U[,which.min(AICfit)],
                  "parameters"=pars[,which.min(AICfit)],
                  "logl"=logl[which.min(AICfit)]))
      
    }
    
  } else if(crit=="BIC"){
    
    BICfit=-2*c(fitN$LogLikelihood,fitLa$LogLikelihood,fitLo$LogLikelihood,fitT$LogLikelihood)+log(n)*c(3,3,3,4)
    
    if(all==T){
      
      return(cbind("margin"=mf,"BIC"=BICfit,"parameters"=t(pars)))
      
    } else {
      
      return(list("margin"=mf[which.min(BICfit)],
                  "U"=U[,which.min(BICfit)],
                  "parameters"=pars[,which.min(BICfit)],
                  "logl"=logl[which.min(BICfit)]))
      
    }
    
    
  } else if(crit=="QQ"){
    
    corfit=c(
      cor(seq(1/(n+1),n/(n+1),length.out=n),sort(U[,1])),
      cor(seq(1/(n+1),n/(n+1),length.out=n),sort(U[,2])),
      cor(seq(1/(n+1),n/(n+1),length.out=n),sort(U[,3])),
      cor(seq(1/(n+1),n/(n+1),length.out=n),sort(U[,4]))
    )
    
    if(all==T){
      
      return(cbind("margin"=mf,"Correlation QQ-plot"=corfit,"parameters"=t(pars)))
      
    } else {
      
      return(list("margin"=mf[which.max(corfit)],
                  "U"=U[,which.max(corfit)],
                  "parameters"=pars[,which.max(corfit)],
                  "logl"=logl[which.max(corfit)]))
      
    }
    
  } else if(crit=="KS"){
    
    suppressWarnings({
      pvalue=c(ks.test(data,pAND,mu=fitN$mu,phi=fitN$phi,alpha=fitN$alpha)$p.value,
               ks.test(data,pALaD,mu=fitLa$mu,phi=fitLa$phi,alpha=fitLa$alpha)$p.value,
               ks.test(data,pALoD,mu=fitLo$mu,phi=fitLo$phi,alpha=fitLo$alpha)$p.value,
               ks.test(data,pATD,mu=fitT$mu,phi=fitT$phi,alpha=fitT$alpha,nu = fitT$nu)$p.value
      )})
    
    if(all==T){
      
      return(cbind("margin"=mf,"KS-test p-value"=pvalue,"parameters"=t(pars)))
      
    } else {
      
      return(list("margin"=mf[which.max(pvalue)],
                  "U"=U[,which.max(pvalue)],
                  "parameters"=pars[,which.max(pvalue)],
                  "logl"=logl[which.max(pvalue)]))
      
    }
  } else {
    
    stop("invalid criterion")
    
  }
  
}  
