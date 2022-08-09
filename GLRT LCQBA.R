### file  for determining the GLRT for symmetry in          ###   
### linear combination of QBA-distributions                 ###
###############################################################


source("~/functions.R")


### LCQBA model parameters ###
##############################

n=1000
d=2
alpha=c(0.5,0.5)
mu=c(0,0)
A=matrix(c(1,0.2,0.2,1),nrow=d,ncol=d)
tpars=c(NA,NA)
basefunc=c("normal","logistic")
seed=sample(x = 1:1000000,size = 1)

X=Xsample(A = A,location = mu,basefunc = basefunc,alpha = alpha,sampsize = n,dims = d,seed = seed,tpars = tpars)
dat=X[[2]]

fitSym=fitLCS(data = dat,basefunc = basefunc,seed = seed,symmetric = T)
fitReg=fitLCS(data = dat,basefunc = basefunc,seed = seed,numstarts = 20,symmetric = F)

Lambda=-2*(fitSym$`log likelihood fit`-fitReg$`log likelihood fit`)
pchisq(q=Lambda,df = d,lower.tail = F)
