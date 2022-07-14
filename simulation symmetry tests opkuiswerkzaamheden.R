source("C:/Users/u0125240/Documents/PhD/code/symmetry/skew normal copula.R")
source("C:/Users/u0125240/Documents/PhD/code/symmetry/symmetry tests.R")
source("C:/Users/u0125240/Documents/PhD/code/lineaire combinatie/Code part 1 Functions.R")
source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitLCS.R")
source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitPcopula.R")
source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitSPcopula.R")
source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitNP.R")

library(copula)
library(ks)
library(sn)
library(nloptr)
library(doParallel)
library(QBAsyDist)
library(simukde)


# function for fitting univariate skew normal distribution to data used in apply
fitUVSN=function(Y,symmetric=T){
  n <- length(Y)
  if(symmetric==T){
    fit <- selm.fit(x = matrix(1,nrow=n,ncol=1),y = Y,family = 'SN',fixed.param = list(alpha=0),selm.control = list(method="MLE"))
    return(fit$param$dp.complete)
  } else {
    fit <- selm.fit(x = matrix(1,nrow=n,ncol=1),y = Y,family = 'SN',selm.control = list(method="MLE"))
    return(fit$param$dp)
  }
}

########################################################################
### setting 1                                                        ###
########################################################################

# load in settings
loadfile = load(paste0("~/PhD/code/symmetry/output/simulation1/modelsettingsAsyDist.Rdata")) # name modelsettings

# dimension and sample size
n <- modelsettings$sampsize
d <- modelsettings$d

# margin parameters
alpha.m <- modelsettings$marginskew
xi.m <- modelsettings$marginlocation
omega.m <- modelsettings$marginscale

# copula parameters (location xi is standard 0)
alpha.c <- modelsettings$copulaskew
rho.c <- modelsettings$copulacor

# number of monte Carlo repetitions
N <- modelsettings$N

# number of starts for optimization in functions
nstart <- 6

# inflation factor bandwidth and grid points non-parametric fit
GP <- modelsettings$GP
IF <- modelsettings$IF

# constants for references tests
a <- modelsettings$a # Henze (2003)
b <- modelsettings$b # Szekely (2001)
c <- modelsettings$c # dai (2018)

indP <- c()
indSP <- c()
for(q in 1:n){
  try({
    load(paste0("~/PhD/code/symmetry/output/simulation1/teststatsP/run",q,".Rdata"))
    indP <- c(indP,q)
    },silent=T)
  try({
    load(paste0("~/PhD/code/symmetry/output/simulation1/teststatsSP/run",q,".Rdata"))
    indSP <- c(indSP,q)
  },silent=T)
}

indPmissing <- c(1:100)[-indP]
indSPmissing <- c(1:100)[-indSP]


### parametric
start.t=Sys.time()

cl <- makeCluster(5)
registerDoParallel(cl)
### main loop for drawing sample from the model and refitting it to the sample
out <- foreach(q=indPmissing,.packages=c('sn','ks','copula','nloptr'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 
                 ### required files:
                 source("C:/Users/u0125240/Documents/PhD/code/symmetry/skew normal copula.R")
                 
                 load(paste0("~/PhD/code/symmetry/output/simulation1/data/run",q,".Rdata")) # dat
                 
                 set.seed(21052022)
                 
                 ### parametric fit assuming margins are fully known, including 
                 ### the parameters
                 try({
                   # regular fit
                   fitReg.P <- fitSNcopula(u = dat$U,symmetric = F,nstart = nstart,random=T)
                   # symmetric fit
                   fitSym.P <- fitSNcopula(u = dat$U,symmetric = T,nstart = nstart,random=T)
                   # test statistic
                   W.P <- -2*(fitSym.P$ll-fitReg.P$ll)
                   # output
                   parametric <- list("fitSym"=fitSym.P,"fitReg"=fitReg.P,
                                      "TS"=W.P)
                   save(parametric,file=paste0("~/PhD/code/symmetry/output/simulation1/teststatsP/run",q,".Rdata"))
                 },silent=TRUE)
               }
stopCluster(cl)
end.t <- Sys.time()
time.takensim <- difftime(end.t,start.t,units = "secs")

### semi-parametric
start.t=Sys.time()

cl <- makeCluster(5)
registerDoParallel(cl)
### main loop for drawing sample from the model and refitting it to the sample
out <- foreach(q=indSPmissing,.packages=c('sn','ks','copula','nloptr'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 
                 ### required files:
                 source("C:/Users/u0125240/Documents/PhD/code/symmetry/skew normal copula.R")
                 
                 load(paste0("~/PhD/code/symmetry/output/simulation1/data/run",q,".Rdata")) # dat
                 
                 set.seed(21052022)
                 
                 ### semi parametric fit and symmetry test using GLRT
                 ### with and without assumed marginal symmetry
                 
                 # generating non-parametric pseudo-observations
                 Pobs.SP <- pobs(x = dat$X)
                 
                 try({
                   # regular fitting
                   fitReg.SP <- fitSNcopula(u = Pobs.SP,symmetric = F,nstart = nstart,random=T)
                   
                   # symmetric fitting with assumed symmetric margins
                   fitSym.SP <- fitSNcopula(u = Pobs.SP,symmetric = T,nstart = nstart,random=T)
                   
                   # test statistic with assumed symmetric margins
                   W.SP <- -2*(fitSym.SP$ll-fitReg.SP$ll)  
                   
                   # output
                   semiparametric <- list("fitSym"=fitSym.SP,"fitReg"=fitReg.SP,
                                          "TS"=W.SP)
                   save(semiparametric,file=paste0("~/PhD/code/symmetry/output/simulation1/teststatsSP/run",q,".Rdata"))
                   
                   
                 },silent=T)
               }
stopCluster(cl)
end.t <- Sys.time()
time.takensim <- difftime(end.t,start.t,units = "secs")

x11()








########################################################################
### setting 2                                                        ###
########################################################################

# load in settings
loadfile = load(paste0("~/PhD/code/symmetry/output/simulation2/modelsettingsAsyDist.Rdata")) # name modelsettings

# dimension and sample size
n <- modelsettings$sampsize
d <- modelsettings$d

# margin parameters
alpha.m <- modelsettings$marginskew
xi.m <- modelsettings$marginlocation
omega.m <- modelsettings$marginscale

# copula parameters (location xi is standard 0)
alpha.c <- modelsettings$copulaskew
rho.c <- modelsettings$copulacor

# number of monte Carlo repetitions
N <- modelsettings$N

# number of starts for optimization in functions
nstart <- 6

# inflation factor bandwidth and grid points non-parametric fit
GP <- modelsettings$GP
IF <- modelsettings$IF

# constants for references tests
a <- modelsettings$a # Henze (2003)
b <- modelsettings$b # Szekely (2001)
c <- modelsettings$c # dai (2018)

indP <- c()
indSP <- c()
for(q in 1:n){
  try({
    load(paste0("~/PhD/code/symmetry/output/simulation2/teststatsP/run",q,".Rdata"))
    indP <- c(indP,q)
  },silent=T)
  # try({
  #   load(paste0("~/PhD/code/symmetry/output/simulation2/teststatsSP/run",q,".Rdata"))
  #   indSP <- c(indSP,q)
  # },silent=T)
}

indPmissing <- c(1:100)[-indP]
# indSPmissing <- c(1:100)[-indSP]


### parametric
start.t=Sys.time()

cl <- makeCluster(5)
registerDoParallel(cl)
### main loop for drawing sample from the model and refitting it to the sample
out <- foreach(q=indPmissing,.packages=c('sn','ks','copula','nloptr'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 
                 ### required files:
                 source("C:/Users/u0125240/Documents/PhD/code/symmetry/skew normal copula.R")
                 
                 load(paste0("~/PhD/code/symmetry/output/simulation2/data/run",q,".Rdata")) # dat
                 
                 set.seed(21052022)
                 
                 ### parametric fit assuming margins are fully known, including 
                 ### the parameters
                 try({
                   # regular fit
                   fitReg.P <- fitSNcopula(u = dat$U,symmetric = F,nstart = nstart,random=T)
                   # symmetric fit
                   fitSym.P <- fitSNcopula(u = dat$U,symmetric = T,nstart = nstart,random=T)
                   # test statistic
                   W.P <- -2*(fitSym.P$ll-fitReg.P$ll)
                   # output
                   parametric <- list("fitSym"=fitSym.P,"fitReg"=fitReg.P,
                                      "TS"=W.P)
                   save(parametric,file=paste0("~/PhD/code/symmetry/output/simulation2/teststatsP/run",q,".Rdata"))
                 },silent=TRUE)
               }
stopCluster(cl)
end.t <- Sys.time()
time.takensim <- difftime(end.t,start.t,units = "secs")

# ### semi-parametric
# start.t=Sys.time()
# 
# cl <- makeCluster(5)
# registerDoParallel(cl)
# ### main loop for drawing sample from the model and refitting it to the sample
# out <- foreach(q=indSPmissing,.packages=c('sn','ks','copula','nloptr'),
#                .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
#                  
#                  ### required files:
#                  source("C:/Users/u0125240/Documents/PhD/code/symmetry/skew normal copula.R")
#                  
#                  load(paste0("~/PhD/code/symmetry/output/simulation2/data/run",q,".Rdata")) # dat
#                  
#                  set.seed(21052022)
#                  
#                  ### semi parametric fit and symmetry test using GLRT
#                  ### with and without assumed marginal symmetry
#                  
#                  # generating non-parametric pseudo-observations
#                  Pobs.SP <- pobs(x = dat$X)
#                  
#                  try({
#                    # regular fitting
#                    fitReg.SP <- fitSNcopula(u = Pobs.SP,symmetric = F,nstart = nstart,random=T)
#                    
#                    # symmetric fitting with assumed symmetric margins
#                    fitSym.SP <- fitSNcopula(u = Pobs.SP,symmetric = T,nstart = nstart,random=T)
#                    
#                    # test statistic with assumed symmetric margins
#                    W.SP <- -2*(fitSym.SP$ll-fitReg.SP$ll)  
#                    
#                    # output
#                    semiparametric <- list("fitSym"=fitSym.SP,"fitReg"=fitReg.SP,
#                                           "TS"=W.SP)
#                    save(semiparametric,file=paste0("~/PhD/code/symmetry/output/simulation2/teststatsSP/run",q,".Rdata"))
#                    
#                    
#                  },silent=T)
#                }
# stopCluster(cl)
# end.t <- Sys.time()
# time.takensim <- difftime(end.t,start.t,units = "secs")
# x11()
























########################################################################
### setting 3                                                        ###
########################################################################

# load in settings
loadfile = load(paste0("~/PhD/code/symmetry/output/simulation3/modelsettingsAsyDist.Rdata")) # name modelsettings

# dimension and sample size
n <- modelsettings$sampsize
d <- modelsettings$d

# margin parameters
alpha.m <- modelsettings$marginskew
xi.m <- modelsettings$marginlocation
omega.m <- modelsettings$marginscale

# copula parameters (location xi is standard 0)
alpha.c <- modelsettings$copulaskew
rho.c <- modelsettings$copulacor

# number of monte Carlo repetitions
N <- modelsettings$N

# number of starts for optimization in functions
nstart <- 6

# inflation factor bandwidth and grid points non-parametric fit
GP <- modelsettings$GP
IF <- modelsettings$IF

# constants for references tests
a <- modelsettings$a # Henze (2003)
b <- modelsettings$b # Szekely (2001)
c <- modelsettings$c # dai (2018)

indP <- c()
indSP <- c()
for(q in 1:n){
  try({
    load(paste0("~/PhD/code/symmetry/output/simulation3/teststatsP/run",q,".Rdata"))
    indP <- c(indP,q)
  },silent=T)
  # try({
  #   load(paste0("~/PhD/code/symmetry/output/simulation3/teststatsSP/run",q,".Rdata"))
  #   indSP <- c(indSP,q)
  # },silent=T)
}

indPmissing <- c(1:100)[-indP]
# indSPmissing <- c(1:100)[-indSP]


### parametric
start.t=Sys.time()

cl <- makeCluster(5)
registerDoParallel(cl)
### main loop for drawing sample from the model and refitting it to the sample
out <- foreach(q=indPmissing,.packages=c('sn','ks','copula','nloptr'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 
                 ### required files:
                 source("C:/Users/u0125240/Documents/PhD/code/symmetry/skew normal copula.R")
                 
                 load(paste0("~/PhD/code/symmetry/output/simulation3/data/run",q,".Rdata")) # dat
                 
                 set.seed(21052022)
                 
                 ### parametric fit assuming margins are fully known, including 
                 ### the parameters
                 try({
                   # regular fit
                   fitReg.P <- fitSNcopula(u = dat$U,symmetric = F,nstart = nstart,random=T)
                   # symmetric fit
                   fitSym.P <- fitSNcopula(u = dat$U,symmetric = T,nstart = nstart,random=T)
                   # test statistic
                   W.P <- -2*(fitSym.P$ll-fitReg.P$ll)
                   # output
                   parametric <- list("fitSym"=fitSym.P,"fitReg"=fitReg.P,
                                      "TS"=W.P)
                   save(parametric,file=paste0("~/PhD/code/symmetry/output/simulation3/teststatsP/run",q,".Rdata"))
                 },silent=TRUE)
               }
stopCluster(cl)
end.t <- Sys.time()
time.takensim <- difftime(end.t,start.t,units = "secs")

# ### semi-parametric
# start.t=Sys.time()
# 
# cl <- makeCluster(5)
# registerDoParallel(cl)
# ### main loop for drawing sample from the model and refitting it to the sample
# out <- foreach(q=indSPmissing,.packages=c('sn','ks','copula','nloptr'),
#                .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
#                  
#                  ### required files:
#                  source("C:/Users/u0125240/Documents/PhD/code/symmetry/skew normal copula.R")
#                  
#                  load(paste0("~/PhD/code/symmetry/output/simulation3/data/run",q,".Rdata")) # dat
#                  
#                  set.seed(21052022)
#                  
#                  ### semi parametric fit and symmetry test using GLRT
#                  ### with and without assumed marginal symmetry
#                  
#                  # generating non-parametric pseudo-observations
#                  Pobs.SP <- pobs(x = dat$X)
#                  
#                  try({
#                    # regular fitting
#                    fitReg.SP <- fitSNcopula(u = Pobs.SP,symmetric = F,nstart = nstart,random=T)
#                    
#                    # symmetric fitting with assumed symmetric margins
#                    fitSym.SP <- fitSNcopula(u = Pobs.SP,symmetric = T,nstart = nstart,random=T)
#                    
#                    # test statistic with assumed symmetric margins
#                    W.SP <- -2*(fitSym.SP$ll-fitReg.SP$ll)  
#                    
#                    # output
#                    semiparametric <- list("fitSym"=fitSym.SP,"fitReg"=fitReg.SP,
#                                           "TS"=W.SP)
#                    save(semiparametric,file=paste0("~/PhD/code/symmetry/output/simulation3/teststatsSP/run",q,".Rdata"))
#                    
#                    
#                  },silent=T)
#                }
# stopCluster(cl)
# end.t <- Sys.time()
# time.takensim <- difftime(end.t,start.t,units = "secs")
# x11()






























########################################################################
### setting 4                                                        ###
########################################################################

source("C:/Users/u0125240/Documents/PhD/code/symmetry/skew normal copula.R")
source("C:/Users/u0125240/Documents/PhD/code/symmetry/symmetry tests.R")
source("C:/Users/u0125240/Documents/PhD/code/lineaire combinatie/Code part 1 Functions.R")
source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitLCS.R")
source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitPcopula.R")
source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitSPcopula.R")
source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitNP.R")

library(copula)
library(ks)
library(sn)
library(nloptr)
library(doParallel)
library(QBAsyDist)
library(simukde)

########################################################################
### setting of the simulation model from the skew-normal copula with ###
### chi-square margins                                               ###
########################################################################

# seed
seed <- 21052022
set.seed(seed)

# number of monte Carlo repetitions
NN <- 400

# load in settings
load(paste0("~/PhD/code/symmetry/output/simulation4/modelsettings.Rdata")) # name modelsettings

# dimension and sample size
n <- modelsettings$sampsize
d <- modelsettings$d

# copula parameters (location xi is standard 0)
alpha.c <- modelsettings$copulaskew
rho.c <- modelsettings$copulacor

# number of monte Carlo repetitions
N <- modelsettings$N

# number of starts for optimization in functions
nstart <- 10

# inflation factor bandwidth and grid points non-parametric fit
GP <- modelsettings$GP
IF <- modelsettings$IF

# constants for references tests
a <- modelsettings$a # Henze (2003)
b <- modelsettings$b # dai (2018)
c <- modelsettings$c # Szekely (2001)

indp <- c()
for(q in 1:N){
  try({
    load(paste0("~/PhD/code/symmetry/output/simulation4/teststatsNP/TSDistRun",q,".Rdata"))
    indp <- c(indp,q)
    
  },silent=T)
}
missind <- c(1:N)[-indp]


# start.t=Sys.time()
# 
# cl <- makeCluster(6)
# registerDoParallel(cl)
# ### main loop for drawing sample from the model and refitting it to the sample
# out <- foreach(q=missind,.packages=c('sn','ks','copula','nloptr'),
#                .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
#                  
#                  ### required files:
#                  source("C:/Users/u0125240/Documents/PhD/code/symmetry/skew normal copula.R")
#                  source("C:/Users/u0125240/Documents/PhD/code/symmetry/symmetry tests.R")
#                  source("C:/Users/u0125240/Documents/PhD/code/lineaire combinatie/Code part 1 Functions.R")
#                  source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitLCS.R")
#                  source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitPcopula.R")
#                  source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitSPcopula.R")
#                  source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitNP.R")
#                  source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitUVSN.R")
#                  
#                  
#                  ### data generation and saving for reference
#                  # creating pseudo-observations
#                  U <- rSNcopula(n = n,alpha = alpha.c,rho = rho.c)
#                  
#                  # creating margins
#                  X <- qchisq(p = U,df = 3)-1
#                  dat <- list("U"=U,"X"=X)
#                  save(dat,file = paste0("~/PhD/code/symmetry/output/simulation4/data/run",q,".Rdata"))
#                  
#                  
#                  
#                  
#                  ### parametric fit and symmetry test using GLRT when margins
#                  ### need to be estimated as well
#                  try({
#                     # regular fitting
#                     fitReg.P.m <- apply(X,MARGIN = 2, FUN = fitUVSN, symmetric=F)
#                     PobsReg <- matrix(NA,nrow=n,ncol=d)
#                     llReg.m <- rep(NA,d)
#                     MPReg <- matrix(NA,nrow=3,ncol=d)
#                     for(i in 1:d){
#                       PobsReg[,i] <- psn(x = X[,i],xi = fitReg.P.m[[i]]$par[1],omega = fitReg.P.m[[i]]$par[2],alpha = fitReg.P.m[[i]]$par[3])
#                       llReg.m[i] <- fitReg.P.m[[i]]$ll
#                       MPReg[,i] <- fitReg.P.m[[i]]$par
#                     }
#                     fitReg.P.c <- fitSNcopula(u = PobsReg,symmetric = F,nstart = nstart)
#                  
#                     # symmetric fitting
#                     fitSym.P.m <- apply(X,MARGIN = 2,FUN = fitUVSN, symmetric=T)
#                     PobsSym <- matrix(NA,nrow=n,ncol=d)
#                     llSym.m <- rep(NA,d)
#                     MPSym <- matrix(NA,nrow=3,ncol=d)
#                     for(i in 1:d){
#                         PobsSym[,i] <- psn(x = X[,i],xi = fitSym.P.m[[i]]$par[1],omega = fitSym.P.m[[i]]$par[2],alpha = fitSym.P.m[[i]]$par[3])
#                     llSym.m[i] <- fitSym.P.m[[i]]$ll
#                     MPSym[,i] <- fitSym.P.m[[i]]$par
#                     }
#                     fitSym.P.c <- fitSNcopula(u = PobsSym,symmetric = T,nstart = nstart)
#                  
#                     # test statistic
#                     W.P_um <- -2*(fitSym.P.c$ll+sum(llSym.m)-fitReg.P.c$ll-sum(llReg.m))
#                  
#                     # output
#                     parametric_um <- list("MPSym"=MPSym,"fitSym"=fitSym.P.c,
#                                           "MPReg"=MPReg,"fitReg"=fitReg.P.c,
#                                           "TS"=W.P_um)
#                     save(parametric_um,file=paste0("~/PhD/code/symmetry/output/simulation4/teststatsP_UM/run",q,".Rdata"))
#                  },silent=T)
#                  
#                  ### parametric fit assuming margins are fully known, including 
#                  ### the parameters
#                  try({
#                    # regular fit
#                    fitReg.P <- fitSNcopula(u = U,symmetric = F,nstart = nstart)
#                    # symmetric fit
#                    fitSym.P <- fitSNcopula(u = U,symmetric = T,nstart = nstart)
#                    # test statistic
#                    W.P <- -2*(fitSym.P$ll-fitReg.P$ll)
#                    # output
#                    parametric <- list("fitSym"=fitSym.P,"fitReg"=fitReg.P,
#                                       "TS"=W.P)
#                    save(parametric,file=paste0("~/PhD/code/symmetry/output/simulation4/teststatsP/run",q,".Rdata"))
#                  },silent=TRUE)
#                  
#                  
#                  ### semi parametric fit and symmetry test using GLRT
#                  ### with and without assumed marginal symmetry
#                  
#                  # generating non-parametric pseudo-observations
#                  Pobs.SP <- pobs(x = X)
#                  
#                  # generating non-parametric symmetrized pseudo-observations 
#                  # with respect to the median 
#                  theta <- apply(X,2,median)
#                  X2 <- sweep(x = -X,MARGIN = 2,STATS = 2*theta,FUN = "+")
#                  Pobs.SPS <- pobs(x = rbind(X,X2))[1:n,]
#                  
#                  try({
#                    # regular fitting
#                    fitReg.SP <- fitSNcopula(u = Pobs.SP,symmetric = F,nstart = nstart,random=T)
#                    
#                    # symmetric fitting with assumed symmetric margins
#                    fitSym.SP <- fitSNcopula(u = Pobs.SP,symmetric = T,nstart = nstart,random=T)
#                    
#                    # test statistic with assumed symmetric margins
#                    W.SP <- -2*(fitSym.SP$ll-fitReg.SP$ll)  
#                    
#                    # output
#                    semiparametric <- list("fitSym"=fitSym.SP,"fitReg"=fitReg.SP,
#                                           "TS"=W.SP)
#                    save(semiparametric,file=paste0("~/PhD/code/symmetry/output/simulation4/teststatsSP/run",q,".Rdata"))
#                  },silent=T)
#                  
#                  
#                  
#                  
#                  ### non-parametric fit and testing using GLRT
#                  ### and reference tests
#                  try({
#                    # regular fitting
#                    fitReg.NP <- fitNP(X = X,symmetric = F,IF = IF,GP = GP)
#                    
#                    # symmetric fitting
#                    fitSym.NP <- fitNP(X = X,symmetric = T,IF = IF,GP = GP)
#                    
#                    # test statistic
#                    W.NP <- -2*(fitSym.NP$ll-fitReg.NP$ll)
#                    
#                    # Dai (2018), Henze (2003) & Szekely (2001) tests
#                    W.DHS <- DHS(X = X,type = "all",a = a,b = b,c = c)
#                    
#                    # Depth based runs test of Dyckerhoff 2015
#                    DRtest <- DBST(X = X,depth = "H")
#                    W.DR <- DRtest$TestValue
#                    
#                    # adding the data augmented with the median reflected data
#                    # for bootstrapping under the hypothesis of symmetry
#                    XAug <- rbind(X,X2)
#                    
#                    # output
#                    nonparametric <- list("bw"= fitSym.NP$bw,"LRTNP"=W.NP,
#                                          "DHS"=W.DHS,"DR"=W.DR,"XAug"=XAug)
#                    save(nonparametric,file=paste0("~/PhD/code/symmetry/output/simulation4/teststatsNP/run",q,".Rdata"))
#                  },silent=T)
#                }
# stopCluster(cl)
# end.t <- Sys.time()
# time.takensim <- difftime(end.t,start.t,units = "secs")

#######################################################################################
### creating the asymptotic distribution under the hypothesis of central symmetry   ###
### for the parametric GLRT, the semi-parametric GLRT using symmetrised margins,    ###
### the non-parametric GLRT and the reference symmetry tests                        ###
#######################################################################################

### non-parametric GLRT and reference tests
###########################################
start.t=Sys.time()

cl <- makeCluster(6)
registerDoParallel(cl)
NPR <- foreach(q=missind,.packages=c('sn','ks','copula','nloptr'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 set.seed(1305222+q)
                 
                 ### required files:
                 source("C:/Users/u0125240/Documents/PhD/code/symmetry/symmetry tests.R")
                 source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitNP.R")
                 
                 # load in the original dataset
                 load(paste0("~/PhD/code/symmetry/output/simulation4/data/run",q,".Rdata"))
                 X <- dat$X
                 
                 # load in the estimates
                 load(paste0("~/PhD/code/symmetry/output/simulation4/teststatsNP/run",q,".Rdata"))
                 bw <- nonparametric$bw 
                 XAug <- nonparametric$XAug
                 
                 # determining the asymptotic distribution using bootstrapping
                 W.NP <- rep(NA,NN)
                 W.DHS <- matrix(NA,nrow=NN,ncol = 3)
                 W.DR <- rep(NA,NN)
                 
                 for(j in 1:NN){
                   # constructing bootstrapped sample
                   samp=jitter(XAug[sample(x = 1:(2*n),size = n,replace = T),],factor = 0.0001)
                   # non-parametric fit and testing using GLRT
                   # and reference tests
                   try({
                     # regular fitting
                     fitReg.NP <- fitNP(X = samp,symmetric = F,IF = IF,GP = GP)
                     # symmetric fitting
                     fitSym.NP <- fitNP(X = samp,symmetric = T,IF = IF,GP = GP)
                     # test statistic
                     W.NP[j] <- -2*(fitSym.NP$ll-fitReg.NP$ll)
                   },silent = T)
                   # Dai (2018), Henze (2003) & Szekely (2001) tests
                   try({
                     W.DHS[j,] <- DHS(X = samp,type = "all",a = a,b = b,c = c)
                   },silent = T)
                   # Depth based runs test of Dyckerhoff 2015
                   try({
                     DRtest <- DBST(X = samp,depth = "H")
                     W.DR[j] <- DRtest$TestValue
                   },silent = T)
                 }         
                 # output
                 TSDist <- list("LRTNP"=W.NP,"DHS"=W.DHS,"DR"=W.DR)          
                 save(TSDist,file=paste0("~/PhD/code/symmetry/output/simulation4/teststatsNP/TSDistRun",q,".Rdata"))
               }
stopCluster(cl)
end.t <- Sys.time()
time.takenstats <- difftime(end.t,start.t,units = "secs")

# ### parametric GLRT with unknown margins
# ########################################
# 
# cl <- makeCluster(6)
# registerDoParallel(cl)
# start.t=Sys.time()
# PR <- foreach(q=missind,.packages=c('sn','ks','copula','nloptr'),
#               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
#                 
#                 set.seed(21052022+q)
#                 
#                 ### required files:
#                 source("C:/Users/u0125240/Documents/PhD/code/symmetry/skew normal copula.R")
#                 source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitUVSN.R")
#                 
#                 
#                 # load in the estimates
#                 load(paste0("~/PhD/code/symmetry/output/simulation4/teststatsP_UM/run",q,".Rdata"))
#                 parCsym <- parametric_um$fitSym$par
#                 parMsym <- parametric_um$MPSym
#                 
#                 
#                 
#                 # determining the asymptotic distribution using Monte Carlo under
#                 # the assumption we don't know the margins 
#                 W.P_um <- rep(NA,NN)
#                 
#                 for(i in 1:NN){
#                   # generate pseudo-observations
#                   U <- rSNcopula(n = n,alpha = rep(0,d),rho = parCsym)
#                   # transform to data scale
#                   X <- matrix(NA,nrow=n,ncol=d)
#                   for(j in 1:d){
#                     X[,j] <- qsn(p = U[,j],xi = parMsym[1,j],omega = parMsym[2,j],alpha = parMsym[3,j])
#                   }
#                   
#                   try({
#                     # regular fitting
#                     fitReg.P.m <- apply(X,MARGIN = 2, FUN = fitUVSN, symmetric=F)
#                     PobsReg <- matrix(NA,nrow=n,ncol=d)
#                     llReg.m <- rep(NA,d)
#                     for(k in 1:d){
#                       PobsReg[,k] <- psn(x = X[,k],xi = fitReg.P.m[[k]]$par[1],omega = fitReg.P.m[[k]]$par[2],alpha = fitReg.P.m[[k]]$par[3])
#                       llReg.m[k] <- fitReg.P.m[[k]]$ll
#                     }
#                     fitReg.P.c <- fitSNcopula(u = PobsReg,symmetric = F,nstart = nstart,random=T)
#                     
#                     # symmetric fitting
#                     fitSym.P.m <- apply(X,MARGIN = 2,FUN = fitUVSN, symmetric=T)
#                     PobsSym <- matrix(NA,nrow=n,ncol=d)
#                     llSym.m <- rep(NA,d)
#                     for(k in 1:d){
#                       PobsSym[,k] <- psn(x = X[,k],xi = fitSym.P.m[[k]]$par[1],omega = fitSym.P.m[[k]]$par[2],alpha = fitSym.P.m[[k]]$par[3])
#                       llSym.m[k] <- fitSym.P.m[[k]]$ll
#                     }
#                     fitSym.P.c <- fitSNcopula(u = PobsSym,symmetric = T,nstart = nstart,random=T)
#                     
#                     # test statistic
#                     W.P_um[i] <- -2*(fitSym.P.c$ll+sum(llSym.m)-fitReg.P.c$ll-sum(llReg.m))
#                   },silent=TRUE)
#                 }         
#                 # output
#                 TSDist <- list("LRTP"=W.P_um)          
#                 save(TSDist,file=paste0("~/PhD/code/symmetry/output/simulation4/teststatsP_UM/TSDistRun",q,".Rdata"))
#               }
# stopCluster(cl)
# end.t <- Sys.time()
# time.takenstats <- difftime(end.t,start.t,units = "secs")
# 
# x11()







#############################
### Redoing the DHS tests ###
#############################

# load in settings
loadfile = load(paste0("~/PhD/code/symmetry/output/simulation1/modelsettingsAsyDist.Rdata")) # name modelsettings

# dimension and sample size
n <- modelsettings$sampsize
d <- modelsettings$d

# margin parameters
alpha.m <- modelsettings$marginskew
xi.m <- modelsettings$marginlocation
omega.m <- modelsettings$marginscale

# copula parameters (location xi is standard 0)
alpha.c <- modelsettings$copulaskew
rho.c <- modelsettings$copulacor

# number of monte Carlo repetitions
N <- modelsettings$N

# number of starts for optimization in functions
nstart <- 6

# inflation factor bandwidth and grid points non-parametric fit
GP <- modelsettings$GP
IF <- modelsettings$IF

# constants for references tests
a <- modelsettings$a # Henze (2003)
b <- modelsettings$b # Szekely (2001)
c <- modelsettings$c # dai (2018)

dir.create("~/PhD/code/symmetry/output/simulation1/teststatsDHS")
start.t=Sys.time()

cl <- makeCluster(5)
registerDoParallel(cl)
### main loop for drawing sample from the model and refitting it to the sample
out <- foreach(q=1:N,.packages=c('sn','ks','copula','nloptr'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 
                 ### required files:
                 source("C:/Users/u0125240/Documents/PhD/code/symmetry/symmetry tests.R")
                 
                 load(paste0("~/PhD/code/symmetry/output/simulation1/data/run",q,".Rdata")) # dat
                 X <- dat$X
                 try({
                   # Dai (2018), Henze (2003) & Szekely (2001) tests
                   W.DHS <- DHS(X = X,type = "all",a = a,b = b,c = c)
                   
                   # output
                   save(W.DHS,file=paste0("~/PhD/code/symmetry/output/simulation1/teststatsDHS/run",q,".Rdata"))
                 },silent=T)
               }
stopCluster(cl)
end.t <- Sys.time()
time.takensim <- difftime(end.t,start.t,units = "secs")


# seed
seed <- 21052022
set.seed(seed)

# number of monte Carlo repetitions
NN <- 400

# load in settings
load(paste0("~/PhD/code/symmetry/output/simulation1/modelsettings.Rdata")) # name modelsettings

# dimension and sample size
n <- modelsettings$sampsize
d <- modelsettings$d

# copula parameters (location xi is standard 0)
alpha.c <- modelsettings$copulaskew
rho.c <- modelsettings$copulacor

# number of monte Carlo repetitions
N <- modelsettings$N

# number of starts for optimization in functions
nstart <- 10

# inflation factor bandwidth and grid points non-parametric fit
GP <- modelsettings$GP
IF <- modelsettings$IF

# constants for references tests
a <- modelsettings$a # Henze (2003)
b <- modelsettings$b # dai (2018)
c <- modelsettings$c # Szekely (2001)

#######################################################################################
### creating the asymptotic distribution under the hypothesis of central symmetry   ###
### for the parametric GLRT, the semi-parametric GLRT using symmetrised margins,    ###
### the non-parametric GLRT and the reference symmetry tests                        ###
#######################################################################################


### non-parametric GLRT and reference tests
###########################################
start.t=Sys.time()

cl <- makeCluster(6)
registerDoParallel(cl)
NPR <- foreach(q=1:N,.packages=c('sn','ks','copula','nloptr'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 set.seed(1305222+q)
                 
                 ### required files:
                 source("C:/Users/u0125240/Documents/PhD/code/symmetry/symmetry tests.R")
                 source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitNP.R")
                 
                 # load in the original dataset
                 load(paste0("~/PhD/code/symmetry/output/simulation1/data/run",q,".Rdata"))
                 X <- dat$X
                 theta <- apply(X,2,median)
                 X2 <- sweep(x = -X,MARGIN = 2,STATS = 2*theta,FUN = "+")
                 XAug <- rbind(X,X2)
                 
                 # determining the asymptotic distribution using bootstrapping
                 W.DHS <- matrix(NA,nrow=NN,ncol = 3)
                 
                 for(j in 1:NN){
                   # constructing bootstrapped sample
                   samp=jitter(XAug[sample(x = 1:(2*n),size = n,replace = T),],factor = 0.0001)
                   # non-parametric fit and testing using GLRT
                   # and reference tests
                   # Dai (2018), Henze (2003) & Szekely (2001) tests
                   try({
                     W.DHS[j,] <- DHS(X = samp,type = "all",a = a,b = b,c = c)
                   },silent = T)
                 }         
                 # output
                 save(W.DHS,file=paste0("~/PhD/code/symmetry/output/simulation1/teststatsDHS/TSDistRun",q,".Rdata"))
               }
stopCluster(cl)
end.t <- Sys.time()
time.takenstats <- difftime(end.t,start.t,units = "secs")

x11()
