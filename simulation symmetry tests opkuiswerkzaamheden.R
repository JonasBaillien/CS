source("C~/functions.R")

library(copula)
library(ks)
library(sn)
library(nloptr)
library(doParallel)
library(QBAsyDist)
library(simukde)

########################################################################
### setting 1                                                        ###
########################################################################

# load in settings
loadfile = load(paste0("~/simulation1/modelsettingsAsyDist.Rdata")) # name modelsettings

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
    load(paste0("~/simulation1/teststatsP/run",q,".Rdata"))
    indP <- c(indP,q)
    },silent=T)
  try({
    load(paste0("~/simulation1/teststatsSP/run",q,".Rdata"))
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
                 source("~/functions.R")
                 
                 load(paste0("~/simulation1/data/run",q,".Rdata")) # data
                 
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
                   save(parametric,file=paste0("~/simulation1/teststatsP/run",q,".Rdata"))
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
                 source("~/functions.R")
                 
                 load(paste0("~/simulation1/data/run",q,".Rdata")) # dat
                 
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
                   save(semiparametric,file=paste0("~/simulation1/teststatsSP/run",q,".Rdata"))
                   
                   
                 },silent=T)
               }
stopCluster(cl)
end.t <- Sys.time()
time.takensim <- difftime(end.t,start.t,units = "secs")




########################################################################
### setting 2                                                        ###
########################################################################

# load in settings
loadfile = load(paste0("~/simulation2/modelsettingsAsyDist.Rdata")) # name modelsettings

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
    load(paste0("~/simulation2/teststatsP/run",q,".Rdata"))
    indP <- c(indP,q)
  },silent=T)
  # try({
  #   load(paste0("~/simulation2/teststatsSP/run",q,".Rdata"))
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
                 source("~/functions.R")
                 
                 load(paste0("~/simulation2/data/run",q,".Rdata")) # dat
                 
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
                   save(parametric,file=paste0("~/simulation2/teststatsP/run",q,".Rdata"))
                 },silent=TRUE)
               }
stopCluster(cl)
end.t <- Sys.time()
time.takensim <- difftime(end.t,start.t,units = "secs")


########################################################################
### setting 3                                                        ###
########################################################################

# load in settings
loadfile = load(paste0("~/simulation3/modelsettingsAsyDist.Rdata")) # name modelsettings

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
    load(paste0("~/simulation3/teststatsP/run",q,".Rdata"))
    indP <- c(indP,q)
  },silent=T)
  # try({
  #   load(paste0("~/simulation3/teststatsSP/run",q,".Rdata"))
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
                 source("~/functions.R")
                 
                 load(paste0("~/simulation3/data/run",q,".Rdata")) # dat
                 
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
                   save(parametric,file=paste0("~/simulation3/teststatsP/run",q,".Rdata"))
                 },silent=TRUE)
               }
stopCluster(cl)
end.t <- Sys.time()
time.takensim <- difftime(end.t,start.t,units = "secs")
