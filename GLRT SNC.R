source("~/functions.R")


### copula model specification ###
##################################
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
