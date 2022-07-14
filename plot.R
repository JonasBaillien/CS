# returns the axes and density to plot in contour(xax,yax,dens)
plotLCS=function(basefunc,alpha,mu,A,tpars,X,n,xlim=NULL,ylim=NULL){
  
  # grid in 1 direction
  if(is.null(xlim)){
    xtemp=seq(min(X[,1])-sqrt(sd(X[,1])),max(X[,1])+sqrt(sd(X[,1])),length.out=n)
  } else {
    xtemp=seq(xlim[1],xlim[2],length.out=n)
  }
  if(is.null(ylim)){
    ytemp=seq(min(X[,2])-sqrt(sd(X[,2])),max(X[,2])+sqrt(sd(X[,2])),length.out=n)
  } else {
    ytemp=seq(ylim[1],ylim[2],length.out=n)
  }
  # square grid
  xx=as.matrix(expand.grid(xtemp,ytemp)) 
  
  fmat=apply(xx,1,densityfLCS,basefunc=basefunc,alpha=alpha,mu=mu,B=solve(A),tpars=tpars)
  fmat=matrix(fmat,nrow=n)
  
  return(list("xax"=xtemp,"yax"=ytemp,"dens"=fmat))
}

densityfLCS=function(point,basefunc,alpha,mu,B,tpars){
  
  d=length(point)
  
  B=matrix(B,nrow=d)
  
  # diag(B)=pmin(10^4,exp(diag(B)))
  
  V=point%*%B
  VV=mu%*%B
  
  dens=abs(det(B))
  
  for(i in 1:d){
    if(basefunc[i]=="laplace"){
      dens=dens*twopiecelaplacedensity(V[i]-VV[i],alpha[i])
    } else if(basefunc[i]=="normal"){
      dens=dens*twopiecenormaldensity(V[i]-VV[i],alpha[i])
    } else if(basefunc[i]=="logistic"){
      dens=dens*twopiecelogisticdensity(V[i]-VV[i],alpha[i])
    } else {
      dens=dens*twopiecestudentdensity(V[i]-VV[i],alpha[i],nu=tpars[i])
    }
  }
  
  return(dens)
}


plotPcopula=function(fitPcop,n,X,xlim=NULL,ylim=NULL){
  
  if(dim(fitPcop$fitC)!=2){
    stop("Plotting only available for dimensions 2")
  } else {
    d=2
  }
  
  # grid in 1 direction
  if(is.null(xlim)){
    xax=seq(min(X[,1])-sqrt(sd(X[,1])),max(X[,1])+sqrt(sd(X[,1])),length.out=n)
  } else {
    xax=seq(xlim[1],xlim[2],length.out=n)
  }
  if(is.null(ylim)){
    yax=seq(min(X[,2])-sqrt(sd(X[,2])),max(X[,2])+sqrt(sd(X[,2])),length.out=n)
  } else {
    yax=seq(ylim[1],ylim[2],length.out=n)
  }
  ax=cbind(xax,yax)
  
  
  densityax=matrix(NA,nrow=n,ncol=d)
  U=matrix(NA,nrow=n,ncol=d)
  
  margfuncs=fitPcop$fitM$margins
  alpha=fitPcop$fitM$parameters[1,]
  mu=fitPcop$fitM$parameters[2,]
  phi=fitPcop$fitM$parameters[3,]
  df=fitPcop$fitM$parameters[4,]
  
  for(i in 1:d){
    if(margfuncs[i]=="normal"){
      U[,i]=pAND(q = ax[,i],mu=mu[i],phi=phi[i],alpha=alpha[i])
      densityax[,i]=dAND(y = ax[,i],mu=mu[i],phi=phi[i],alpha=alpha[i])
    } else if(margfuncs[i]=="logistic"){
      U[,i]=pALoD(q = ax[,i],mu=mu[i],phi=phi[i],alpha=alpha[i])
      densityax[,i]=dALoD(y = ax[,i],mu=mu[i],phi=phi[i],alpha=alpha[i])
    } else if(margfuncs[i]=="laplace"){
      U[,i]=pALaD(q = ax[,i],mu=mu[i],phi=phi[i],alpha=alpha[i])
      densityax[,i]=dALaD(y = ax[,i],mu=mu[i],phi=phi[i],alpha=alpha[i])
    } else {
      U[,i]=pATD(q = ax[,i],mu=mu[i],phi=phi[i],alpha=alpha[i],nu=df[i])
      densityax[,i]=dATD(y = ax[,i],mu=mu[i],phi=phi[i],alpha=alpha[i],nu=df[i])
    }
  }
  
  Ugrid=as.matrix(expand.grid(U[,1],U[,2]))
  densitycopula=matrix(dCopula(Ugrid,fitPcop$fitC@copula),nrow=n) # copula density on the [0,1]^2 grid
  
  dens=matrix(NA,nrow=n,ncol=n) # density on the data scale
  for(i in 1:n){
    for(j in 1:n){
      dens[i,j]=densityax[i,1]*densityax[j,2]*densitycopula[i,j]
    }
  }
  
  return(list("xax"=xax,"yax"=yax,"dens"=dens))
  
}

