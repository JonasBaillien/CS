library(copula)
library(QBAsyDist)
source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitPcopula.R")
source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitLCS.R")
source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitSPcopula.R")
source("C:/Users/u0125240/Documents/PhD/code/symmetry/skew normal copula.R")
source("C:/Users/u0125240/Documents/PhD/code/symmetry/plot.R")

n=100
d=3
alpha=c(0.5,0.5,0.2)
mu=c(0,-2,2)
phi=c(0.5,1,2)
sigma=c(0.5,0.2,-0.4)

cop=normalCopula(param = sigma,dim = d,dispstr = "un")
U=rCopula(copula = cop,n = n)

X=matrix(NA,n,d)
X[,1]=qAND(beta = U[,1],mu = mu[1],phi = phi[1],alpha = alpha[1])
X[,2]=qALoD(beta = U[,2],mu = mu[2],phi = phi[2],alpha = alpha[2])
X[,3]=qALaD(beta = U[,3],mu = mu[3],phi = phi[3],alpha = alpha[3])


fitcop=normalCopula(dim = d,dispstr = "un")
margfuncs=c("normal","logistic","laplace")

f1=fitCopulaQBAM(X = X,cop = fitcop,margfunc = NULL,symmetric = F)
f2=fitCopulaQBAM(X = X,cop = fitcop,margfunc = NULL,symmetric = T)
f3=fitCopulaQBAM(X = X,cop = fitcop,margfunc = margfuncs,symmetric = F)
f4=fitCopulaQBAM(X = X,cop = fitcop,margfunc = margfuncs,symmetric = T)

library(DAAG)
ais=DAAG::ais[,c(6,9)]
f1=fitCopulaQBAM(X = ais,cop = normalCopula(dim = 2,dispstr = "un"),symmetric = F)
f2=fitCopulaQBAM(X = ais,cop = normalCopula(dim = 2,dispstr = "un"),symmetric = T)
f3=fitSPcopula(X = ais,cop = normalCopula(dim=2,dispstr = "un" ),theta.method = "median",symmetric = F)
f4=fitSPcopula(X = ais,cop = normalCopula(dim=2,dispstr = "un" ),theta.method = "median",symmetric = T)
f5=fitLCS(data = ais,basefunc = c("logistic","normal"),symmetric = F)
f6=fitLCS(data = ais,basefunc = c("logistic","normal"),symmetric = T)





# load data below in the matrix X
f3=fitCopulaQBAM(X =  X,cop = tCopula(dim = 2,dispstr = "un"),symmetric = F)
f4=fitCopulaQBAM(X =  X,cop = tCopula(dim = 2,dispstr = "un"),symmetric = T)

p3=plotPcopula(fitPcop = f3,n = 500,X = X,xlim = c(50,300),ylim=c(20,50))
p4=plotPcopula(fitPcop = f4,n = 500,X = X,xlim = c(50,300),ylim=c(20,50))
x11()
par(mar=c(5.1, 5.1, 4.1, 2.1))
contour(p3$xax,p3$yax,p3$dens,
        #xlab=expression(X[1]),ylab=expression(X[2]),
        #xlab="Weight",ylab="Waist",
        #xlim=c(-0.07,0.07),ylim=c(-0.07,0.07),
        cex.lab=1.5,cex.axis=1.5,lwd=2,cex=2)
contour(p4$xax,p4$yax,p4$dens,col=2,lwd=2,add=T)




basefunc=c("t","logistic")
f1=fitLCS(data = X,basefunc = basefunc,symmetric = F,seed = 1245)
f2=fitLCS(data = X,basefunc = basefunc,symmetric = T,seed = 1245)

p1=plotLCS(basefunc = basefunc,alpha = as.numeric(f1$`fitted parameters`[1:2]),
        mu = as.numeric(f1$`fitted parameters`[3:4]),
        A = matrix(as.numeric(f1$`fitted parameters`[5:8]),nrow = 2),
        tpars=as.numeric(f1$`fitted parameters`[9]),X = X,n = 200,xlim = c(50,300),ylim=c(20,50))    

p2=plotLCS(basefunc = basefunc,alpha = as.numeric(f2$`fitted parameters`[1:2]),
        mu = as.numeric(f2$`fitted parameters`[3:4]),
        A = matrix(as.numeric(f2$`fitted parameters`[5:8]),nrow = 2),
        tpars = as.numeric(f2$`fitted parameters`[9]),X = X,n = 200,xlim = c(50,300),ylim=c(20,50))
x11()
par(mar=c(5.1, 5.1, 4.1, 2.1))
contour(p1$xax,p1$yax,p1$dens,
        #xlab=expression(X[1]),ylab=expression(X[2]),
        xlab="Weight",ylab="Waist",
        cex.lab=1.5,cex.axis=1.5,lwd=2,cex=2)
contour(p2$xax,p2$yax,p2$dens,col=2,lwd=2,add=T,levels=seq(2e-05,0.0012,by=2e-004))



x11()
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(X,cex.lab=1.5,cex.axis=1.5,lwd=2,cex=1,pch=19,xlim=c(-0.07,0.07),ylim=c(-0.07,0.07),
     xlab="DAX",ylab="CAC40")


### empirical density estimation ###
####################################
library(ks)
n=1000
xmin = c(50,20)
xmax = c(300,50)
fittest=kde(x = X,gridsize = c(n,n),xmin = xmin,xmax = xmax)
bw=2*fittest$H
fit=kde(x = X,H = bw,gridsize = c(n,n),xmin = xmin,xmax = xmax)
theta.ind=which(fit$estimate == max(fit$estimate),arr.ind = T)
theta=c(fit$eval.points[[1]][theta.ind[1]],fit$eval.points[[2]][theta.ind[2]])
fitm=kde(x = as.matrix(sweep(x = -X,MARGIN = 2,STATS = 2*theta,FUN = "+")),H = bw,gridsize = c(n,n),xmin = xmin,xmax = xmax)
xgrid=fit$eval.points[[1]]
ygrid=fit$eval.points[[2]]

x11()
contour(xgrid,ygrid,fit$estimate,cex.lab=1.5,cex.axis=1.5,lwd=2,
        xlab="Weigth",ylab="Waist",levels=seq(0.0002,0.0014,by=0.0002))
contour(xgrid,ygrid,0.5*(fit$estimate+fitm$estimate),add=T,col=2,cex.lab=1.5,cex.axis=1.5,lwd=2,levels=seq(0.0002,0.0014,by=0.0002))


n=1000
xmin = c(-0.07,-0.07)
xmax = c(0.07,0.07)
fittest=kde(x = X,gridsize = c(n,n),xmin = xmin,xmax = xmax)
bw=2*fittest$H
fit=kde(x = X,H = bw,gridsize = c(n,n),xmin = xmin,xmax = xmax)
theta.ind=which(fit$estimate == max(fit$estimate),arr.ind = T)
theta=c(fit$eval.points[[1]][theta.ind[1]],fit$eval.points[[2]][theta.ind[2]])
fitm=kde(x = as.matrix(sweep(x = -X,MARGIN = 2,STATS = 2*theta,FUN = "+")),H = bw,gridsize = c(n,n),xmin = xmin,xmax = xmax)
xgrid=fit$eval.points[[1]]
ygrid=fit$eval.points[[2]]

x11()
contour(xgrid,ygrid,fit$estimate,cex.lab=1.5,cex.axis=1.5,lwd=2,
        xlab="DAX",ylab="CAC40",levels=seq(100,800,by=100))
contour(xgrid,ygrid,0.5*(fit$estimate+fitm$estimate),add=T,col=2,cex.lab=1.5,
        cex.axis=1.5,lwd=2,levels=seq(100,800,by=100))

#### Data examples ###
######################
wineR <- read.csv("~/PhD/code/datasets/winequality-red.csv", sep=";")
table(wineR$quality)
x11()
pairs(wineR[wineR$quality==6,c(2,9,10)])
# te onderzoeken opties
# quality 6: Volatile acidity en pH
# quality 6: Sulphates en pH
X=as.matrix(wineR[wineR$quality==6,c(2,10)])




################

indi <- read.csv("~/PhD/code/datasets/indian_liver_patient.csv")
x11()
pairs(indi[,8:10])
# total proteins en albumin en albumin globulin ratio
X=na.omit(indi[,c(8,10)])
######################

wineW <- read.csv("~/PhD/code/datasets/winequality-white.csv", sep=";")
table(wineW$quality)
x11()
pairs(wineW[wineW$quality==4,c(2,3,8,9)])
x11()
pairs(wineW[wineW$quality==5,])
x11()
pairs(wineW[wineW$quality==6,c(1,9:11)])
# meerdere opties


####################

# in lessR library
body <- read.csv("~/PhD/code/datasets/dataBodyMeas.csv")
x11()
pairs(body[,-2])
X=body[,c(3,5)]

####################

data(pbc, package = "randomForestSRC")
x11()
pairs(pbc[,c(11,12,16)]) # also 12 with values <2.6 removed also 16 log transformed
x11()
pairs(pbc[,c(10:18)])

#######################

library(quantmod)
getSymbols(c('^GDAXI','^FCHI','^N225','AMZN','NFLX','AMD'),from = "2010-01-01", to = "2018-12-31")
CAC40=as.numeric(weeklyReturn(FCHI,type="log"))[-470]
DAX=as.numeric(weeklyReturn(GDAXI,type="log"))
Nikkei225=as.numeric(weeklyReturn(N225,type="log"))
Amazon=as.numeric(weeklyReturn(AMZN,type="log"))
Netflix=as.numeric(weeklyReturn(NFLX,type="log"))
AMD=as.numeric(weeklyReturn(AMD,type="log"))

pairs(cbind(CAC40,DAX,Nikkei225,Amazon,Netflix,AMD))

X=cbind(DAX,CAC40)



