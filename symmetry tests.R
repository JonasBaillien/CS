########################
### Test of symmetry ###
########################
# D and S need to be swapped
library(MTS)
DHS=function(X,type="all",a=2,b=1,c=1){
  # X: nxd matrix containing the data
  # type: "H" (Henze 2003), "S" (Szekely 2001), "D" (Dai 2018) or "all" (default)
  # a: constant of Henze 2003 to govern power in (0,infty]. Defaults to 2
  # b: constant of Szekely 2001 to govern power in (0,2]. Defaults to 1
  # c: constant of Dai 2018 to govern power in (0,2]. Defaults to 1
  
  # dimensions of data
  n=nrow(X)
  d=ncol(X)
  
  # permutations
  U=sample(x = c(-1,1),size = n,replace = T,prob = c(0.5,0.5))
  
  # calculating Z
  Av=colMeans(X)
  Sn=cov(X)
  Y = sweep(x = X,MARGIN = 2,STATS = Av,FUN = "-")
  Z = sweep(x = sweep(x = sweep(x = X,MARGIN = 2,STATS = Av,FUN = "-"),MARGIN = 1,STATS = U,FUN = "*")%*%msqrt(solve(Sn))$mtxsqrt,MARGIN = 1,STATS = sample(x = c(-1,1),size = n,replace = T),FUN = "*")
  AvZ=colMeans(Z)
  if(type=="H"){
    # test statistic of Henze (21)
    YH=0
    for(i in 1:n){
      for(j in 1:n){
        YH=YH+2+sum(AvZ^2)/(a)-(1+(Z[i,]-Z[j,])%*%AvZ/(2*a))^2*exp(-sum((Z[i,]-Z[j,])^2)/(4*a)) -
          (1+(Z[i,]+Z[j,])%*%AvZ/(2*a))^2*exp(-sum((Z[i,]+Z[j,])^2)/(4*a))
      }
    }
    YH=YH*pi^(d/2)/(2*n*a^(d/2))
    out=c(YH)
    names(out)=c("H")
    
  } else if(type=="S"){
    # test statistic of Szekely (22)
    YS=0
    for(i in 1:n){
      for(j in 1:n){
        YS=YS+exp(-sum((Y[i,]-Y[j,])^2)^(b/2)) -
          exp(-sum((Y[i,]+Y[j,])^2)^(b/2))
      }
    }
    YS=1/(n^2)*YS
    out=c(YS)
    names(out)=c("S")
    
  } else if(type=="D"){
    # test statistic of Dai (23)
    YD=0
    for(i in 1:n){
      for(j in 1:n){
        YD=YD-sum((Y[i,]-Y[j,])^2)^(c/2) +
          sum((Y[i,]+Y[j,])^2)^(c/2)
      }
    }
    YD=1/(n^2)*YD
    out=c(YD)
    names(out)=c("D")
    
  } else if(type=="all"){
    # all three of the above
    
    # test statistic of Henze (21)
    YH=0
    for(i in 1:n){
      for(j in 1:n){
        YH=YH+2+sum(Av^2)/(a)-(1+(Z[i,]-Z[j,])%*%Av/(2*a))^2*exp(-sum((Z[i,]-Z[j,])^2)/(4*a)) -
          (1+(Z[i,]+Z[j,])%*%Av/(2*a))^2*exp(-sum((Z[i,]+Z[j,])^2)/(4*a))
      }
    }
    YH=YH*pi^(d/2)/(2*n*a^(d/2))
    
    # test statistic of Szekely (23)
    YS=0
    for(i in 1:n){
      for(j in 1:n){
        YS=YS+exp(-sum((Y[i,]-Y[j,])^2)^(b/2)) -
          exp(-sum((Y[i,]+Y[j,])^2)^(b/2))
      }
    }
    YS=1/(n^2)*YS

    # test statistic of Dai (22)
    YD=0
    for(i in 1:n){
      for(j in 1:n){
        YD=YD-sum((Y[i,]-Y[j,])^2)^(c/2) +
          sum((Y[i,]+Y[j,])^2)^(c/2)
      }
    }
    YD=1/(n^2)*YD
    out=c(YH,YS,YD)
    names(out)=c("H","S","D")
  }
  return(out) 
}

# example  
# library(mvtnorm)
# X=rmvnorm(n = 100,mu = rep(0,2),sigma = diag(2))
# DHS(X = X,type = "all",a = 4,b = 1.5,c = 1)

#############################
### Depth based runs test ###
#############################
library(geometry)
library(mrfDepth)
library(mvtnorm)
DBST=function(X,depth="SP"){
  # X: nxd matrix with observations in the row
  # depth: depth function to be used; "H"=halspace, "S"=simplicial,
  #        "P"=projection, "DP"=directional projection,
  #        "SP"= skew adjusted projection depth
  
  n=nrow(X)
  d=ncol(X)
  
  if(d==2 & depth=="S"){
    error("Simplicial depth is only implemented for d=2")
  }
  
  switch(depth,
         "H"={
             # estimating the center
             theta=hdepthmedian(x = X)$median
             # calculating depth of data in augmented data
             Xbar=rbind(X,sweep(x = -X,MARGIN = 2,STATS = 2*theta,FUN = "+"))
             # depth
             depth=hdepth(x = Xbar,z = X)$depthZ
         },
         "S"={
           theta=hdepthmedian(x = X)$median
           Xbar=rbind(X,sweep(x = -X,MARGIN = 2,STATS = 2*theta,FUN = "+"))
           depth=sdepth(x = Xbar,z = X)$depthZ
         },
         "P"={
           theta=projmedian(x = X)$gravity
           Xbar=rbind(X,sweep(x = -X,MARGIN = 2,STATS = 2*theta,FUN = "+"))
           depth=projdepth(x = Xbar,z = X)$depthZ
         },
         "DP"={
           theta=dprojmedian(x = X)$gravity
           Xbar=rbind(X,sweep(x = -X,MARGIN = 2,STATS = 2*theta,FUN = "+"))
           depth=dprojdepth(x = Xbar,z = X)$depthZ
         },
         "SP"={
           theta=sprojmedian(x = X)$gravity
           Xbar=rbind(X,sweep(x = -X,MARGIN = 2,STATS = 2*theta,FUN = "+"))
           depth=sprojdepth(x = Xbar,z = X)$depthZ
         })
  
  # sort data according by descending depth
  Z=cbind(X,depth)
  Z=Z[order(Z[,d+1],decreasing = T),-(d+1)]
  TestValue=1
  for(i in (d+1):n){
    ch=convhulln(p = Z[i:(i-d),],output.options = T)
    TestValue=TestValue+1*inhulln(ch = ch,p = matrix(theta,nrow=1,ncol=d)) 
  }
  
  if(d==2){
    pvalue=pnorm(q = (4*TestValue-n-2)/sqrt(n*11/3),mean = 0,sd = 1,lower.tail = T)
    
    return(list("TestValue"=TestValue,"P-value"=pvalue,"Median"=theta))
  } else {
    return(list("TestValue"=TestValue,"Median"=theta))
  }
  
}

# examples
# X=rmvnorm(n = 200,mu = rep(0,3),sigma = diag(3))
# DBST(X = X,depth = "H")
# TV=rep(NA,1000)
# Med=matrix(NA,nrow=1000,ncol=3)
# Sigma=matrix(c(20,0.5,-1,0.5,1,0.9,-1,0.9,20),nrow=3,ncol=3)
# for(i in 1:1000){  
#   X=rmvt(n = 100,delta = rep(4,3),sigma = Sigma,df = 2)
#   #X=rmvnorm(n = 100,mu = rep(4,3),sigma = Sigma)
#   #X=rmsn(n = 100,xi = rep(4,3),Omega = Sigma,alpha = c(10,10,-10))
#   out=DBST(X = X,depth="SP")
#   out
#   TV[i]=out$TestValue
#   Med[i,]=out$Median
#   print(i)
# }


