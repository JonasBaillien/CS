###############################################################################################
### This is for simulation setting 3, for others, change the parameters and the directories ###
###############################################################################################

source("functions.R")

library(copula)
library(ks)
library(sn)
library(nloptr)
library(doParallel)
library(QBAsyDist)
library(simukde)

########################################################################
### setting of the simulation model from the skew-normal copula with ###
### skew-normal margins                                              ###
########################################################################

# seed
seed <- 25422
set.seed(seed)

# dimension and sample size
n <- 100
d <- 3

# margin parameters
alpha.m <- c(0,0,0) # skew setting 1 & 2
# alpha.m <- c(5,1,-3) # setting 3
xi.m <- c(4,3,-2) # location
omega.m <- c(2,0.5,1) # scale
 
# copula parameters (location xi is standard 0)
alpha.c <- c(0,0,0) # skew setting 1
# alpha.c <- c(3,-2,2) # setting 2 & 3
rho.c <- c(0.2,0.1,-0.3) # correlations, length = d*(d-1)/2, all settings

# number of repetitions
# number of monte Carlo repetitions
N <- 100


# create folder to store output in
dir.create(paste0("simulation3"))


#######################################################################
### creating samples and getting test statistics for all used tests ###
#######################################################################

# number of starts for optimization in functions
nstart <- 4

# inflation factor bandwidth and grid points non-parametric fit
GP <- 100
IF <- 2

# constants for references tests
a <- 2 # Henze (2003)
b <- 1 # dai (2018)
c <- 1 # Szekely (2001)

modelsettings <- list("seed"=seed,"d"=d,"sampsize"=n,"marginskew"=alpha.m,
                      "marginlocation"=xi.m,"marginscale"=omega.m,
                      "copulaskew"=alpha.c,"copulacor"=rho.c,"N"=N,
                      "nstart"=nstart,"GP"=GP,"IF"=IF,"a"=a,"b"=b,"c"=c)
save(modelsettings,file=paste0("simulation3/modelsettings.Rdata"))

# creating directories to save output
dir.create("simulation3/teststatsP")
dir.create("teststatsSP")
dir.create("teststatsSPS")
dir.create("teststatsNP")
dir.create("simulation3/data")


start.t=Sys.time()

cl <- makeCluster(5)
registerDoParallel(cl)
### main loop for drawing sample from the model and refitting it to the sample
out <- foreach(q=1:N,.packages=c('sn','ks','copula','nloptr'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 
                 ### required files:
                 source("functions.R")
                 
                 
                 ### data generation and saving for reference
                 # creating pseudo-observations
                 U <- rSNcopula(n = n,alpha = alpha.c,rho = rho.c)
                 
                 # creating margins 
                 X <- matrix(NA,nrow=n,ncol=d)
                 for(i in 1:d){
                   X[,i] <- qsn(p = U[,i],xi = xi.m[i],omega = omega.m[i],alpha = alpha.m[i])
                 }
                 dat <- list("U"=U,"X"=X)
                 save(dat,file = paste0("simulation3/data/run",q,".Rdata"))


                 
                 
                 ### parametric fit assuming margins are fully known, including
                 ### the parameters
                 try({
                 # regular fit
                 fitReg.P <- fitSNcopula(u = U,symmetric = F,nstart = nstart)
                 # symmetric fit
                 fitSym.P <- fitSNcopula(u = U,symmetric = T,nstart = nstart)
                 # test statistic
                 W.P <- -2*(fitSym.P$ll-fitReg.P$ll)
                 # output
                 parametric <- list("fitSym"=fitSym.P,"fitReg"=fitReg.P,
                                    "TS"=W.P)
                 save(parametric,file=paste0("simulation3/teststatsP/run",q,".Rdata"))
                 },silent=TRUE)
                 
                 
                 ### non-parametric fit and testing using GLRT
                 ### and reference tests
                 theta <- NULL
                 try({
                 # regular fitting
                 fitReg.NP <- fitNP(X = X,symmetric = F,IF = IF,GP = GP)
                 
                 # symmetric fitting
                 fitSym.NP <- fitNP(X = X,symmetric = T,IF = IF,GP = GP)
                 
                 # test statistic
                 W.NP <- -2*(fitSym.NP$ll-fitReg.NP$ll)
                 
                 # Dai (2018), Henze (2003) & Szekely (2001) tests
                 W.DHS <- DHS(X = X,type = "all",a = a,b = b,c = c)

                 # Depth based runs test of Dyckerhoff 2015
                 DRtest <- DBST(X = X,depth = "H")
                 W.DR <- DRtest$TestValue
                 
                 # adding the data augmented with the median reflected data
                 # for bootstrapping under the hypothesis of symmetry
                 # generating augmented data
                 theta <- fitSym.NP$mode
                 X2 <- sweep(x = -X,MARGIN = 2,STATS = 2*theta,FUN = "+")
                 XAug <- rbind(X,X2)
                 
                 # output
                 nonparametric <- list("bw"= fitSym.NP$bw,"LRTNP"=W.NP,"XAug"=XAug)#,"DHS"=W.DHS,"DR"=W.DR)
                 save(nonparametric,file=paste0("simulation3/teststatsNP/runNPGLRT",q,".Rdata"))
                 },silent=T)
                 
                 
                 ### semi parametric fit and symmetry test using GLRT
                 ### with and without assumed marginal symmetry
                 
                 # generating non-parametric pseudo-observations
                 Pobs.SP <- pobs(x = X)
                 try({
                 if(is.null(theta)){
                 theta <- apply(X,2,median)
                 X2 <- sweep(x = -X,MARGIN = 2,STATS = 2*theta,FUN = "+")
                 Pobs.SPS <- pobs(x = rbind(X,X2))[1:n,]
                 }
                 # regular fitting
                 fitReg.SP <- fitSNcopula(u = Pobs.SP,symmetric = F,nstart = nstart,random=T)

                 # symmetric fitting with assumed symmetric margins
                 fitSym.SP <- fitSNcopula(u = Pobs.SP,symmetric = T,nstart = nstart,random=T)

                 # test statistic with assumed symmetric margins
                 W.SP <- -2*(fitSym.SP$ll-fitReg.SP$ll)

                 # output
                 semiparametric <- list("fitSym"=fitSym.SP,"fitReg"=fitReg.SP,
                                        "TS"=W.SP)
                 save(semiparametric,file=paste0("simulation3/teststatsSP/run",q,".Rdata"))


                 # symmetric fitting without assumed symmetric margins
                 fitSym.SPS <- fitSNcopula(u = Pobs.SPS,symmetric = T,nstart = nstart,random=T)

                 # test statistic without  assumed symmetric margins
                 W.SPS <- -2*(fitSym.SPS$ll-fitReg.SP$ll)

                 # output
                 semiparametric2 <- list("fitSym"=fitSym.SPS,"fitReg"=fitReg.SP,
                                    "TS"=W.SPS)
                 save(semiparametric2,file=paste0("simulation3/teststatsSPS/run",q,".Rdata"))
                 },silent=T)
}
stopCluster(cl)
end.t <- Sys.time()
time.takensim <- difftime(end.t,start.t,units = "secs")




#######################################################################################
### creating the asymptotic distribution under the hypothesis of central symmetry   ###
### for the semi-parametric GLRT using symmetrised margins, the non-parametric GLRT ###
### and the reference symmetry tests                                                ###
#######################################################################################

# number of monte Carlo repetitions
NN <- 400

# load in settings
load(paste0("simulation3/modelsettings.Rdata")) # name modelsettings

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
nstart <- modelsettings$nstart

# inflation factor bandwidth and grid points non-parametric fit
GP <- modelsettings$GP
IF <- modelsettings$IF

# constants for references tests
a <- modelsettings$a # Henze (2003)
b <- modelsettings$b # Szekely (2001)
c <- modelsettings$c # dai (2018)


# save new settings
modelsettings["NN"] <- NN
save(modelsettings,file = paste0("simulation3/modelsettingsAsyDist.Rdata"))


### non-parametric GLRT and reference tests
###########################################
start.t=Sys.time()

cl <- makeCluster(5)
registerDoParallel(cl)
NPR <- foreach(q=1:N,.packages=c('sn','ks','copula','nloptr'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 
                 set.seed(27422+q)
                 
                 ### required files:
                 source("functions.R")
                 
                 # load in the original dataset
                 load(paste0("simulation3/data/run",q,".Rdata"))
                 X <- dat$X
                 
                 # load in the estimates
                 load(paste0("simulation3/teststatsNP/run",q,".Rdata"))
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
                 TSDist <- list("LRTNP"=W.NP)#,"DHS"=W.DHS,"DR"=W.DR)          
                 save(TSDist,file=paste0("simulation3/teststatsNP/TSDistRun",q,".Rdata"))
                 
              }
stopCluster(cl)
end.t <- Sys.time()
time.takenstats <- difftime(end.t,start.t,units = "secs")


###############################
### calculation of p-values ###
###############################


### parametric testing under the assumption the margins are known & 
### semi-parametric with assumed symmetric margins
pvalue.P <- rep(NA,N)
pvalue.SP <- rep(NA,N)
for(q in 1:N){
  # parametric
  # load in the estimates
  load(paste0("simulation3/teststatsP/run",q,".Rdata"))
  # p-value
  pvalue.P[q] <- pchisq(q = parametric$TS,df = d,lower.tail = FALSE)
  
  load(paste0("simulation3/teststatsSP/run",q,".Rdata"))
  pvalue.SP[q]<- pchisq(q = semiparametric$TS,df = d,lower.tail = FALSE )
}


### non-parametric testing and reference tests
pvalue.NP <- rep(NA,N)
pvalue.H <- rep(NA,N)
pvalue.S <- rep(NA,N)
pvalue.D <- rep(NA,N)
pvalue.R <- rep(NA,N)
for(q in 1:N){
  
  # load in the estimates
  load(paste0("simulation3/teststatsNP/TSDistRun",q,".Rdata"))
  load(paste0("simulation3/teststatsNP/run",q,".Rdata"))
  
  # p-values
  pvalue.NP[q] <- sum(nonparametric$LRTNP<TSDist$LRTNP)/NN
  pvalue.H[q] <- sum(nonparametric$DHS[1]<TSDist$DHS[,1])/NN
  pvalue.S[q] <- sum(nonparametric$DHS[2]<TSDist$DHS[,2])/NN
  pvalue.D[q] <- sum(nonparametric$DHS[3]<TSDist$DHS[,3])/NN
  pvalue.R[q] <- sum(nonparametric$DR<TSDist$DR)/NN
  
}
