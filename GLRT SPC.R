source("~/PhD/code/symmetry/fitSPcopula.R")
library(QBAsyDist)
# model specification
d=2
sigma=c(0.4)#,0.5,-0.5)
cop=ellipCopula(family = "normal",param = sigma,dim = d,dispstr = "un")

# generating pseudo observations
n=1000
datU=rCopula(n = n,copula = cop)

# transforming to data plane
dat=matrix(NA,nrow=n,ncol=d)
dat[,1]=qAND(beta = datU[,1],mu = 2,phi = 2,alpha = 0.35)
dat[,2]=qALoD(beta = datU[,2],mu = -2,phi = 0.5,alpha = 0.5)
#dat[,3]=qALoD(beta = datU[,3],mu = 4,phi = 0.1,alpha = 0.5)

# fit semi-parametric model to data under assumption of symmetry
FitSym=fitSPcopula(X = dat,cop = cop,theta.method = "median",symmetric = T)
LLSym=sum(dCopula(u = pobs(dat),copula = FitSym@copula,log = T))

# fit semi-parametric model to data
FitReg=fitSPcopula(X=dat,cop = cop,symmetric = F)
LLReg=sum(dCopula(u = pobs(dat),copula = FitReg@copula,log = T))

# GLRT statistic
W=-2*(LLSym-LLReg)

N=1000
WMC=rep(NA,N)
set.seed(12)
for(i in 1:N){
  U=rCopula(n = n,copula = FitSym@copula)
  U2=sweep(x = -U,MARGIN = 2,STATS = rep(1,d),"+")
  Us=pobs(rbind(U,U2))[1:n,]
  
  fitS=fitCopula(normalCopula(dim = d,dispstr = "un"),data = Us)
  LLS=sum(dCopula(u = U,copula = fitS@copula,log = T))
  fitA=fitCopula(normalCopula(dim = d,dispstr = "un"),data = U)
  WMC[i]=-2*(LLS-fitA@loglik)
  print(i)
}
hist(WMC)
abline(v=W,col=2)
P=sum(WMC>W)/N
