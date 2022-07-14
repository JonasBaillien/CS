library(copula)
library(QBAsyDist)
library(ks)
source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitLCS.R")
source("C:/Users/u0125240/Documents/PhD/code/symmetry/skew normal copula.R")
source("C:/Users/u0125240/Documents/PhD/code/symmetry/plot.R")
source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitNP.R")
source("C:/Users/u0125240/Documents/PhD/code/symmetry/QBAmarginfit.R")
source("C:/Users/u0125240/Documents/PhD/code/symmetry/symmetry tests.R")

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



### Body measurement data ###
#############################

body <- read.csv("~/PhD/code/datasets/dataBodyMeas.csv")
X1=as.matrix(body[,c(3,5)])

x11()
par(mar=c(5.1, 5.1, 4.1, 2.1),mai=c(1,1,0.1,0.1))
plot(X1,cex.lab=1.5,cex.axis=1.5,lwd=2,
        xlab="Weigth",ylab="Waist",pch=19,
        ylim=c(15,70),xlim=c(40,380))


### skew-normal copula with QBA-margins ###
###########################################
seed=3158
set.seed(seed)

n1=nrow(X1)
d1=ncol(X1)

# determining margins
margins1R.1=marginfit(data = as.matrix(X1[,1]),crit = "AIC",symmetric=F) # t
margins1R.2=marginfit(data = as.matrix(X1[,2]),crit = "AIC",symmetric=F) # logistic

margins1S.1=marginfit(data = as.matrix(X1[,1]),crit = "AIC",symmetric=T) # t
margins1S.2=marginfit(data = as.matrix(X1[,2]),crit = "AIC",symmetric=T) # t

# pseudo-observations
UP1.R = cbind(margins1R.1$U,margins1R.2$U)
UP1.S = cbind(margins1S.1$U,margins1S.2$U)

# fits
fit1P.R = fitSNcopula(u = UP1.R,symmetric = F,nstart = 5,random = F) 
fit1P.S = fitSNcopula(u = UP1.S,symmetric = T,nstart = 5,random = T) 

### plotting
# grid in 1 direction
N=800
xax=seq(00,400,length.out=N)
yax=seq(00,80,length.out=N)

# regular
densityax=matrix(NA,nrow=N,ncol=d1)
U=matrix(NA,nrow=N,ncol=d1)
# symmetric
densityaxS=matrix(NA,nrow=N,ncol=d1)
US=matrix(NA,nrow=N,ncol=d1)

# regular
U[,1]=pATD(q = xax,mu=margins1R.1$parameters[2],
           phi=margins1R.1$parameters[3],alpha=margins1R.1$parameters[1],
           nu=margins1R.1$parameters[4])
densityax[,1]=dATD(y = xax,mu=margins1R.1$parameters[2],
                   phi=margins1R.1$parameters[3],alpha=margins1R.1$parameters[1],
                   nu=margins1R.1$parameters[4])
U[,2]=pALoD(q = yax,mu=margins1R.2$parameters[2],
           phi=margins1R.2$parameters[3],alpha=margins1R.2$parameters[1])
densityax[,2]=dALoD(y = yax,mu=margins1R.2$parameters[2],
                   phi=margins1R.2$parameters[3],alpha=margins1R.2$parameters[1])
Ugrid=as.matrix(expand.grid(U[,1],U[,2]))
densitycopula=matrix(dSNcopula(u = Ugrid,alpha = fit1P.R$par[1:2],rho = fit1P.R$par[3]),nrow = N,ncol=N) # copula density on the [0,1]^2 grid

dens=matrix(NA,nrow=N,ncol=N) # density on the data scale
for(i in 1:N){
  for(j in 1:N){
    dens[i,j]=densityax[i,1]*densityax[j,2]*densitycopula[i,j]
  }
}

# symmetric
US[,1]=pATD(q = xax,mu=margins1S.1$parameters[2],
           phi=margins1S.1$parameters[3],alpha=margins1S.1$parameters[1],
           nu=margins1S.1$parameters[4])
densityaxS[,1]=dATD(y = xax,mu=margins1S.1$parameters[2],
                   phi=margins1S.1$parameters[3],alpha=margins1S.1$parameters[1],
                   nu=margins1S.1$parameters[4])
US[,2]=pATD(q = yax,mu=margins1S.2$parameters[2],
            phi=margins1S.2$parameters[3],alpha=margins1S.2$parameters[1],
            nu=margins1S.2$parameters[4])
densityaxS[,2]=dATD(y = yax,mu=margins1S.2$parameters[2],
                   phi=margins1S.2$parameters[3],alpha=margins1S.2$parameters[1],
                   nu=margins1S.2$parameters[4])
UgridS=as.matrix(expand.grid(US[,1],US[,2]))
densitycopulaS=matrix(dSNcopula(u = UgridS,alpha = rep(0,2),rho = fit1P.S$par[1]),nrow = N,ncol=N) # copula density on the [0,1]^2 grid
densS=matrix(NA,nrow=N,ncol=N) # density on the data scale
for(i in 1:N){
  for(j in 1:N){
    densS[i,j]=densityaxS[i,1]*densityaxS[j,2]*densitycopulaS[i,j]
  }
}

x11()
par(mar=c(5.1, 5.1, 4.1, 2.1),mai=c(1,1,0.1,0.1))
contour(xax,yax,dens,
        xlab="Weight",ylab="Waist",
        cex.lab=1.5,cex.axis=1.5,lwd=2,cex=2,
        levels=seq(2e-05,0.0012,by=2e-004),
        ylim=c(15,70),xlim=c(40,380))
contour(xax,yax,densS,col=2,lwd=2,add=T,
        levels=seq(2e-05,0.0012,by=2e-004))




### linear combination ###
##########################
basefunc1=c("t","logistic")
basefunc2=c("t","t")
N=1000
f1=fitLCS(data = X1,basefunc = basefunc1,symmetric = F,seed = 1245)
f2=fitLCS(data = X1,basefunc = basefunc2,symmetric = T,seed = 1245)

p1=plotLCS(basefunc = basefunc1,alpha = as.numeric(f1$`fitted parameters`[1:2]),
           mu = as.numeric(f1$`fitted parameters`[3:4]),
           A = matrix(as.numeric(f1$`fitted parameters`[5:8]),nrow = 2),
           tpars=as.numeric(f1$`fitted parameters`[9]),X = X1,n = N,xlim = c(0,400),ylim=c(0,80))    

p2=plotLCS(basefunc = basefunc2,alpha = as.numeric(f2$`fitted parameters`[1:2]),
           mu = as.numeric(f2$`fitted parameters`[3:4]),
           A = matrix(as.numeric(f2$`fitted parameters`[5:8]),nrow = 2),
           tpars = as.numeric(f2$`fitted parameters`[9:10]),X = X1,n = N,xlim = c(0,400),ylim=c(0,80))
x11()
par(mar=c(5.1, 5.1, 4.1, 2.1))
contour(p1$xax,p1$yax,p1$dens,
        levels=seq(2e-05,0.0012,by=2e-004),
        xlab="Weight",ylab="Waist",
        cex.lab=1.5,cex.axis=1.5,lwd=2,cex=2,
        ylim=c(15,70),xlim=c(40,380))
contour(p2$xax,p2$yax,p2$dens,col=2,lwd=2,levels=seq(2e-05,0.0012,by=2e-004),add=T)

### symmetry test
W1.LC <- -2*(f2$`log likelihood fit`-f1$`log likelihood fit`)
pvalue1.LC <- pchisq(q = W1.LC,df=2,lower.tail = F)


### non-parametric fit ###
##########################
n=1000
xmin = c(0,0)
xmax = c(400,80)
fittest=kde(x = X1,gridsize = c(n,n),xmin = xmin,xmax = xmax)
bw=2*fittest$H
fit=kde(x = X1,H = bw,gridsize = c(n,n),xmin = xmin,xmax = xmax)
theta.ind=which(fit$estimate == max(fit$estimate),arr.ind = T)
theta=c(fit$eval.points[[1]][theta.ind[1]],fit$eval.points[[2]][theta.ind[2]])
fitm=kde(x = as.matrix(sweep(x = -X1,MARGIN = 2,STATS = 2*theta,FUN = "+")),H = bw,gridsize = c(n,n),xmin = xmin,xmax = xmax)
xgrid=fit$eval.points[[1]]
ygrid=fit$eval.points[[2]]

x11()
contour(xgrid,ygrid,fit$estimate,cex.lab=1.5,cex.axis=1.5,lwd=2,
        xlab="Weigth",ylab="Waist",levels=seq(2e-05,0.0012,by=2e-004),
        ylim=c(15,70),xlim=c(40,380))
contour(xgrid,ygrid,0.5*(fit$estimate+fitm$estimate),add=T,col=2,cex.lab=1.5,cex.axis=1.5,lwd=2,levels=seq(0.0002,0.0014,by=0.0002))


### symmetry tests for copula and non-parametric ###
####################################################
body <- read.csv("~/PhD/code/datasets/dataBodyMeas.csv")
X1=as.matrix(body[,c(3,5)])
seed=3158
set.seed(seed)

n1=nrow(X1)
d1=ncol(X1)
### parametric fitting
# determining margins
margins1R.1=marginfit(data = as.matrix(X1[,1]),crit = "AIC",symmetric=F) # t
margins1R.2=marginfit(data = as.matrix(X1[,2]),crit = "AIC",symmetric=F) # logistic

margins1S.1=marginfit(data = as.matrix(X1[,1]),crit = "AIC",symmetric=T) # t
margins1S.2=marginfit(data = as.matrix(X1[,2]),crit = "AIC",symmetric=T) # logistic

# pseudo-observations
UP1.R = cbind(margins1R.1$U,margins1R.2$U)
UP1.S = cbind(margins1S.1$U,margins1S.2$U)

# fits
fit1P.R = fitSNcopula(u = UP1.R,symmetric = F,nstart = 10,random = T) 
fit1P.S = fitSNcopula(u = UP1.S,symmetric = T,nstart = 10,random = T) 


### semi-parametric fitting

# symmetrized sample w.r.t. the median
theta1 <- apply(X1,2,median)
X1S <- sweep(x = -X1,MARGIN = 2,STATS = 2*theta1,FUN = "+")
X1AUG <- rbind(X1,X1S)

# pseudo-observations
USP1 = pobs(X1)

# regular
fit1SP.R = fitSNcopula(u = USP1,symmetric = F,nstart = 10,random = T)
# symmetric on regular pseudo-observations
fit1SP.S = fitSNcopula(u = USP1,symmetric = T,nstart = 10,random = T)

### non-parametric fitting
fit1NP.R = fitNP(X = X1,symmetric = F,IF = 2,GP = 100,ub = c(550,75),lb = c(50,15))
fit1NP.S = fitNP(X = X1,symmetric = T,IF = 2,GP = 100,ub = c(550,75),lb = c(-200,0))

### testing for symmetry
W1.P <- -2*(fit1P.S$ll+margins1S.1$logl+margins1S.2$logl-fit1P.R$ll-margins1R.1$logl-margins1R.2$logl)
W1.SP <- -2*(fit1SP.S$ll-fit1SP.R$ll) # direct with chi-square
W1.NP <- -2*(fit1NP.S$ll-fit1NP.R$ll)
W1.DHS <- DHS(X = X1,type = "all",a = 2,b = 1,c = 1)
W1.DR <- DBST(X = X1+runif(n1*d1,-0.001,0.001),depth = "H") # N(n+2,11/3n)

DW1.P <- rep(NA,400)
DW1.NP <- rep(NA,400)
DW1.D <- rep(NA,400)
DW1.H <- rep(NA,400)
DW1.S  <- rep(NA,400)
DW1.DR <- rep(NA,400)

for(i in 1:400){
  load(paste0("~/PhD/code/symmetry/output/dataexamples/bodymeas/AsymptDistBodyMeasRun",i,".Rdata"))
  DW1.D[i] <- out$MCW1.DHS[2]
  DW1.S[i] <- out$MCW1.DHS[3]
  DW1.H[i] <- out$MCW1.DHS[1]
  DW1.DR[i] <- out$MCW1.DR$TestValue
  load(paste0("~/PhD/code/symmetry/output/dataexamples/bodymeas/AsymptDistPBodyMeasRun",i,".Rdata"))
  DW1.P[i] <- out$MCW1.P
  load(paste0("~/PhD/code/symmetry/output/dataexamples/bodymeas/AsymptDistNPRun",i,".Rdata"))
  DW1.NP[i] <- out$MCW1.NP
}

pvalue1.P <- sum(DW1.P>W1.P)/400
pvalue1.SP <- pchisq(q = W1.SP,df = 2,lower.tail=F)
pvalue1.NP <- sum(DW1.NP>W1.NP)/400
pvalue1.S <- sum(DW1.S>W1.DHS[3])/400
pvalue1.H <- sum(DW1.H>W1.DHS[1])/400
pvalue1.D <- sum(DW1.D<W1.DHS[2])/400
pvalue1.DR <- W1.DR$`P-value`
# sum(DW1.DR<W1.DR$TestValue)/400













##################
### Stock data ###
##################

library(quantmod)
getSymbols(c('^GDAXI','^FCHI'),from = "2010-01-01", to = "2018-12-31")
CAC40=as.numeric(weeklyReturn(FCHI,type="log"))[-470]
DAX=as.numeric(weeklyReturn(GDAXI,type="log"))
X2=as.matrix(cbind(DAX,CAC40))

x11()
par(mar=c(5.1, 5.1, 4.1, 2.1),mai=c(1,1,0.1,0.1))
plot(X2,xlab="DAX",ylab="CAC40",pch=19,
     cex.lab=1.5,cex.axis=1.5,lwd=2,cex=1,
     xlim=c(-0.07,0.07),ylim=c(-0.07,0.07))

### skew-normal copula with QBA-margins ###
###########################################
seed=3158
set.seed(seed)

n2=nrow(X2)
d2=ncol(X2)

### parametric fitting
# determining margins
margins2R.1=marginfit(data = as.matrix(X2[,1]),crit = "AIC",symmetric=F) # logistic
margins2R.2=marginfit(data = as.matrix(X2[,2]),crit = "AIC",symmetric=F) # logistic

margins2S.1=marginfit(data = as.matrix(X2[,1]),crit = "AIC",symmetric=T) # logistic
margins2S.2=marginfit(data = as.matrix(X2[,2]),crit = "AIC",symmetric=T) # logistic

# pseudo-observations
UP2.R = cbind(margins2R.1$U,margins2R.2$U)
UP2.S = cbind(margins2S.1$U,margins2S.2$U)

# fits
fit2P.R = fitSNcopula(u = UP2.R,symmetric = F,nstart = 5,random = T) 
fit2P.S = fitSNcopula(u = UP2.S,symmetric = T,nstart = 5,random = T) 


### plotting
# grid in 1 direction
N=400
xax=seq(-0.07,0.07,length.out=N)
yax=seq(-0.07,0.07,length.out=N)

# regular
densityax=matrix(NA,nrow=N,ncol=d2)
U=matrix(NA,nrow=N,ncol=d2)
# symmetric
densityaxS=matrix(NA,nrow=N,ncol=d2)
US=matrix(NA,nrow=N,ncol=d2)

# regular
U[,1]=pALoD(q = xax,mu=margins2R.1$parameters[2],
           phi=margins2R.1$parameters[3],alpha=margins2R.1$parameters[1])
densityax[,1]=dALoD(y = xax,mu=margins2R.1$parameters[2],
                   phi=margins2R.1$parameters[3],alpha=margins2R.1$parameters[1])
U[,2]=pALoD(q = yax,mu=margins2R.2$parameters[2],
            phi=margins2R.2$parameters[3],alpha=margins2R.2$parameters[1])
densityax[,2]=dALoD(y = yax,mu=margins2R.2$parameters[2],
                    phi=margins2R.2$parameters[3],alpha=margins2R.2$parameters[1])
Ugrid=as.matrix(expand.grid(U[,1],U[,2]))
densitycopula=matrix(dSNcopula(u = Ugrid,alpha = fit2P.R$par[1:2],rho = fit2P.R$par[3]),nrow = N,ncol=N) # copula density on the [0,1]^2 grid

dens=matrix(NA,nrow=N,ncol=N) # density on the data scale
for(i in 1:N){
  for(j in 1:N){
    dens[i,j]=densityax[i,1]*densityax[j,2]*densitycopula[i,j]
  }
}

# symmetric
US[,1]=pALoD(q = xax,mu=margins2S.1$parameters[2],
            phi=margins2S.1$parameters[3],alpha=margins2S.1$parameters[1])
densityaxS[,1]=dALoD(y = xax,mu=margins2S.1$parameters[2],
                    phi=margins2S.1$parameters[3],alpha=margins2S.1$parameters[1])
US[,2]=pALoD(q = yax,mu=margins2S.2$parameters[2],
             phi=margins2S.2$parameters[3],alpha=margins2S.2$parameters[1])
densityaxS[,2]=dALoD(y = yax,mu=margins2S.2$parameters[2],
                     phi=margins2S.2$parameters[3],alpha=margins2S.2$parameters[1])
UgridS=as.matrix(expand.grid(US[,1],US[,2]))
densitycopulaS=matrix(dSNcopula(u = UgridS,alpha = rep(0,2),rho = fit2P.S$par[1]),nrow = N,ncol=N) # copula density on the [0,1]^2 grid
densS=matrix(NA,nrow=N,ncol=N) # density on the data scale
for(i in 1:N){
  for(j in 1:N){
    densS[i,j]=densityaxS[i,1]*densityaxS[j,2]*densitycopulaS[i,j]
  }
}

x11()
par(mar=c(5.1, 5.1, 4.1, 2.1),mai=c(1,1,0.1,0.1))
contour(xax,yax,dens,
        xlab="DAX",ylab="CAC40",
        cex.lab=1.5,cex.axis=1.5,lwd=2,cex=2,
        levels=seq(100,800,by=100))
contour(xax,yax,densS,col=2,lwd=2,add=T,levels=seq(100,800,by=100))


### linear combination ###
##########################
basefunc=c("logistic","logistic")
N=400
f3=fitLCS(data = X2,basefunc = basefunc,symmetric = F,seed = 1245,numstarts = 5)
f4=fitLCS(data = X2,basefunc = basefunc,symmetric = T,seed = 1245,numstarts = 5)

p3=plotLCS(basefunc = basefunc,alpha = as.numeric(f3$`fitted parameters`[1:2]),
           mu = as.numeric(f3$`fitted parameters`[3:4]),
           A = matrix(as.numeric(f3$`fitted parameters`[5:8]),nrow = 2),
           X = X2,n = N,xlim = c(-0.07,0.07),ylim=c(-0.07,0.07))    

p4=plotLCS(basefunc = basefunc,alpha = as.numeric(f4$`fitted parameters`[1:2]),
           mu = as.numeric(f4$`fitted parameters`[3:4]),
           A = matrix(as.numeric(f4$`fitted parameters`[5:8]),nrow = 2),
           X = X2,n = N,xlim = c(-0.07,0.07),ylim=c(-0.07,0.07))
x11()
par(mar=c(5.1, 5.1, 4.1, 2.1))
contour(p3$xax,p3$yax,p3$dens,
        xlab="DAX",ylab="CAC40",
        cex.lab=1.5,cex.axis=1.5,lwd=2,cex=2,
        levels=seq(100,800,by=100))
contour(p4$xax,p4$yax,p4$dens,col=2,lwd=2,add=T,levels=seq(100,800,by=100))

### symmetry test
W2.LC <- -2*(f4$`log likelihood fit`-f3$`log likelihood fit`)
pvalue2.LC <- pchisq(q = W2.LC,df=2,lower.tail = F)

### non-parametric fit ###
##########################
N=400
xmin = c(-0.07,-0.07)
xmax = c(0.07,0.07)
fittest=kde(x = X2,gridsize = c(N,N),xmin = xmin,xmax = xmax)
bw=2*fittest$H
fit=kde(x = X,H = bw,gridsize = c(N,N),xmin = xmin,xmax = xmax)
theta.ind=which(fit$estimate == max(fit$estimate),arr.ind = T)
theta=c(fit$eval.points[[1]][theta.ind[1]],fit$eval.points[[2]][theta.ind[2]])
fitm=kde(x = as.matrix(sweep(x = -X,MARGIN = 2,STATS = 2*theta,FUN = "+")),H = bw,gridsize = c(N,N),xmin = xmin,xmax = xmax)
xgrid=fit$eval.points[[1]]
ygrid=fit$eval.points[[2]]

x11()
contour(xgrid,ygrid,fit$estimate,cex.lab=1.5,cex.axis=1.5,lwd=2,
        xlab="DAX",ylab="CAC40",levels=seq(100,800,by=100))
contour(xgrid,ygrid,0.5*(fit$estimate+fitm$estimate),add=T,col=2,cex.lab=1.5,
        cex.axis=1.5,lwd=2,levels=seq(100,800,by=100))

### symmetry tests for copula and non-parametric ###
####################################################

### parametric fitting
# determining margins
margins2R.1=marginfit(data = as.matrix(X2[,1]),crit = "AIC",symmetric=F) # logistic
margins2R.2=marginfit(data = as.matrix(X2[,2]),crit = "AIC",symmetric=F) # logistic

margins2S.1=marginfit(data = as.matrix(X2[,1]),crit = "AIC",symmetric=T) # logistic
margins2S.2=marginfit(data = as.matrix(X2[,2]),crit = "AIC",symmetric=T) # logistic

# pseudo-observations
UP2.R = cbind(margins2R.1$U,margins2R.2$U)
UP2.S = cbind(margins2S.1$U,margins2S.2$U)

# fits
fit2P.R = fitSNcopula(u = UP2.R,symmetric = F,nstart = 10,random = T) 
fit2P.S = fitSNcopula(u = UP2.S,symmetric = T,nstart = 10,random = T) 


### semi-parametric fitting

# pseudo-observations
USP2 = pobs(X2)

# regular
fit2SP.R = fitSNcopula(u = USP2,symmetric = F,nstart = 10,random = T)
# symmetric on regular pseudo-observations
fit2SP.S = fitSNcopula(u = USP2,symmetric = T,nstart = 10,random = T)

### non-parametric fitting
fit2NP.R = fitNP(X = X2,symmetric = F,IF = 2,GP = 100,ub = c(0.16,0.135),lb = c(-0.15,-0.13))
fit2NP.S = fitNP(X = X2,symmetric = T,IF = 2,GP = 100,ub = c(0.16,0.135),lb = c(-0.15,-0.13))

### testing for symmetry
W2.P <- -2*(fit2P.S$ll+margins2S.1$logl+margins2S.2$logl-fit2P.R$ll-margins2R.1$logl-margins2R.2$logl)
W2.SP <- -2*(fit2SP.S$ll-fit2SP.R$ll) # direct with chi-square
W2.NP <- -2*(fit2NP.S$ll-fit2NP.R$ll)
W2.DHS <- DHS(X = X2,type = "all",a = 2,b = 1,c = 1)
W2.DR <- DBST(X = X2+runif(n1*d1,-0.000001,0.0000001),depth = "H") # N(n+2,11/3n)

DW2.P <- rep(NA,400)
DW2.NP <- rep(NA,400)
DW2.D <- rep(NA,400)
DW2.H <- rep(NA,400)
DW2.S  <- rep(NA,400)
DW2.DR <- rep(NA,400)

for(i in 1:400){
  load(paste0("~/PhD/code/symmetry/output/dataexamples/stock/AsymptDistStockRun",i,".Rdata"))
  DW2.NP[i] <- out$MCW2.NP
  DW2.D[i] <- out$MCW2.DHS[2]
  DW2.S[i] <- out$MCW2.DHS[3]
  DW2.H[i] <- out$MCW2.DHS[1]
  DW2.DR[i] <- out$MCW2.DR$TestValue
  load(paste0("~/PhD/code/symmetry/output/dataexamples/stock/AsymptDistPStockRun",i,".Rdata"))
  DW2.P[i] <- out$MCW2.P
}

pvalue2.P <- sum(DW2.P>W2.P)/400
pvalue2.SP <- pchisq(q = W2.SP,df = 2,lower.tail=F)
pvalue2.NP <- sum(DW2.NP>W2.NP)/400
pvalue2.S <- sum(DW2.S>W2.DHS[3])/400
pvalue2.H <- sum(DW2.H>W2.DHS[1])/400
pvalue2.D <- sum(DW2.D<W2.DHS[2])/400
pvalue2.DR <- W2.DR$`P-value`
# sum(DW2.DR<W2.DR$TestValue)/400
