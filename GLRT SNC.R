source("~/PhD/code/symmetry/skew normal copula.R")

library(sn)
library(nloptr)
library(doParallel)

# copula model specification
d=2
rho=c(0.2)#,0.5,-0.5)
alpha=c(1,2)
nstart=3
N=100


# generating pseudo observations
n=1000
datU=rSNcopula(n = n,alpha = alpha,rho = rho)

# fit skew-normal copula to data under assumption of symmetry
FitSym=fitSNcopula(u = datU,nstart = nstart,symmetric = T)
LLSym=FitSym$ll

# fit skew-normal copula to data
FitReg=fitSNcopula(u = datU,nstart = nstart,symmetric = F)
LLReg=FitReg$ll

# GLRT statistic
W=-2*(LLSym-LLReg)
qchisq(p = 0.95,df = 1)
# in Yoshiba 2018 they also use the chi square for the LRT when checking
# for significant skewness parameters



# ### simulating data and fitting the model to it
# cores=detectCores()
# cl <- makeCluster(6)
# registerDoParallel(cl)
# 
# ### create folder to store output in
# dir.create(paste0("~/PhD/code/symmetry/output/testSNcop"))
# 
# ### settings
# settings=list("seed"=seed,"d"=d,"sampsize"=n,"nstart"=nstart,"alpha"=alpha,"rho"=rho,"reps"=N)
# save(settings,file=paste0("~/PhD/code/symmetry/output/testSNcop_settings.Rdata"))
# 
# 
# ### main loop for drawing sample from the model and refitting it to the sample
# result=foreach(i=1:N,.packages=c('sn'),
#                .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
#                  
#                  ### required files:
#                  source("~/PhD/code/symmetry/skew normal copula.R")
#                  
#                  ### data generation and fitting
#                  U=rSNcopula(n = n,alpha = rep(0,d),rho = FitSym$par)
#                  
#                  fitS=fitSNcopula(u = U,symmetric = T,nstart = nstart)
#                  fitA=fitSNcopula(u = U,symmetric = F,nstart = nstart)
#                  WMC=-2*(fitS$ll-fitA$ll)
#                  
#                  ### output
#                  fit=list("fitSym"=fitS,"fitReg"=fitA,"Wstat"=WMC)
#                  
#                  filename=paste0("~/PhD/code/symmetry/output/testSNcop/run",i,".Rdata")
#                  save(fit,file=filename)
#                }
# 
# stopCluster(cl)
# 
# # Bundling the output
# dataset=list()
# ind=c()
# for (i in 1:N){
#   try({
#     load(paste0("~/PhD/code/symmetry/output/testSNcop/run",i,".Rdata"))
#     dataset[[i]] <- fit
#     ind=c(ind,i)
#   },silent=T)
# }
# save(dataset,file=paste0("~/PhD/code/symmetry/output/testSNcop.Rdata"))
# wcm=rep(NA,N)
# for(i in ind){
#   wcm[i]=dataset[[i]]$Wstat
# }
# wcm=na.omit(wcm)
# S=length(ind)
# 
# hist(wcm)
# abline(v=W,col=2)
# Pvalue=sum(wcm>W)/S

library(sn)
library(doParallel)

# copula model specification
d=2
n=1000
rho=c(0.3)
alpha=c(0,0)
nstart=4
N=2000
seed=18
set.seed(seed)
### simulating data and fitting the model to it
cores=detectCores()
cl <- makeCluster(6)
registerDoParallel(cl)

### create folder to store output in
dir.create(paste0("~/PhD/code/symmetry/output/testSNcopchisq2"))

### settings
settings=list("seed"=seed,"d"=d,"sampsize"=n,"nstart"=nstart,"rho"=rho,"reps"=N)
save(settings,file=paste0("~/PhD/code/symmetry/output/testSNcopchisq2_settings.Rdata"))


### main loop for drawing sample from the model and refitting it to the sample
result=foreach(i=1:N,.packages=c('sn'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {

                 ### required files:
                 source("~/PhD/code/symmetry/skew normal copula.R")

                 ### data generation and fitting
                 U=rSNcopula(n = n,alpha = alpha,rho = rho)

                 fitS=fitSNcopula(u = U,symmetric = T,nstart = nstart)
                 fitA=fitSNcopula(u = U,symmetric = F,nstart = nstart)
                 WMC=-2*(fitS$ll-fitA$ll)
                 
                 reject=1*(WMC>qchisq(p = 0.95,df = d))
                 ### output
                 fit=list("fitSym"=fitS,"fitReg"=fitA,"Wstat"=WMC,"reject"=reject)

                 filename=paste0("~/PhD/code/symmetry/output/testSNcopchisq2/run",i,".Rdata")
                 save(fit,file=filename)
               }

stopCluster(cl)

# Bundling the output
dataset=list()
ind=c()
for (i in 1:N){
  try({
    load(paste0("~/PhD/code/symmetry/output/testSNcopchisq2/run",i,".Rdata"))
    dataset[[i]] <- fit
    ind=c(ind,i)
  },silent=T)
}
save(dataset,file=paste0("~/PhD/code/symmetry/output/testSNcopchisS2.Rdata"))
rmc=rep(NA,N)
for(i in ind){
  rmc[i]=dataset[[i]]$reject
}
table(rmc)
