### load in required functions ###
##################################
source("~/functions.R")


#### Data examples ###
######################

### Example 1: Body measeruments data
# loading data (from lessR package)
library(lessR)
data("dataBodyMeas")
body <- dataBodyMeas

# plotting data
x11()
pairs(body[,-2])
X1=as.matrix(body[,c(3,5)])

### fitting a skew-normal copula with QBA-margins ###
#####################################################
seed=3158
set.seed(seed)

n1=nrow(X1)
d1=ncol(X1)

### parametric fitting

## determining marginals
# unconstrained
margins1R.1=marginfit(data = as.matrix(X1[,1]),crit = "AIC",symmetric=F) # t
margins1R.2=marginfit(data = as.matrix(X1[,2]),crit = "AIC",symmetric=F) # logistic

# symmetric
margins1S.1=marginfit(data = as.matrix(X1[,1]),crit = "AIC",symmetric=T) # t
margins1S.2=marginfit(data = as.matrix(X1[,2]),crit = "AIC",symmetric=T) # logistic

## pseudo-observations
# unconstrained
UP1.R = cbind(margins1R.1$U,margins1R.2$U)
# symmetric
UP1.S = cbind(margins1S.1$U,margins1S.2$U)


## fits
# unconstrained
fit1P.R = fitSNcopula(u = UP1.R,symmetric = F,nstart = 10,random = T) 
# symmetric
fit1P.S = fitSNcopula(u = UP1.S,symmetric = T,nstart = 10,random = T) 



### semi-parametric fitting

## symmetrized sample w.r.t. the median
theta1 <- apply(X1,2,median)
X1S <- sweep(x = -X1,MARGIN = 2,STATS = 2*theta1,FUN = "+")
X1AUG <- rbind(X1,X1S)

## pseudo-observations
# empirical original data
USP1 = pobs(X1)
# empirical in augmented data
USPS1 = pobs(X1AUG)[1:n1,]

## fitting SN-copula
# regular
fit1SP.R = fitSNcopula(u = USP1,symmetric = F,nstart = 10,random = T)
# symmetric on regular pseudo-observations
fit1SP.S = fitSNcopula(u = USP1,symmetric = T,nstart = 10,random = T)
# symmetric on symmetrized pseudo-observations
fit1SPS.S = fitSNcopula(u = USPS1,symmetric = T,nstart = 10,random = T)

### non-parametric fitting
# unconstrained
fit1NP.R = fitNP(X = X1,symmetric = F,IF = 2,GP = 100,ub = c(550,75),lb = c(50,15))
# centrally symmetric
fit1NP.S = fitNP(X = X1,symmetric = T,IF = 2,GP = 100,ub = c(550,75),lb = c(-200,0))



### testing for symmetry
W1.P <- -2*(fit1P.S$ll-fit1P.R$ll)
W1.SP <- -2*(fit1SP.S$ll-fit1SP.R$ll) # direct with chi-square
W1.SPS <- -2*(fit1SPS.S$ll-fit1SP.R$ll)
W1.NP <- -2*(fit1NP.S$ll-fit1NP.R$ll)
W1.DHS <- DHS(X = X1,type = "all",a = 2,b = 1,c = 1)
W1.DR <- DBST(X = X1+runif(n1*d1,-0.001,0.001),depth = "H")

### Asymptotic distribution of test statistics
dir.create("~/dataexamples")

testvalues=list("W1.P"=W1.P,"W1.SP"=W1.SP,"W1.SPS"=W1.SPS,"W1.NP"=W1.NP,"W1.DHS"=W1.DHS,"W1.DR"=W1.DR)
save(testvalues,file=paste0("~/dataexamples/BodyMeasTestValues.Rdata"))

# parametric tests in parallel
cl <- makeCluster(5)
registerDoParallel(cl)
### main loop for drawing sample from the model and refitting it to the sample
out <- foreach(q=1:400,.packages=c('sn','ks','copula','nloptr'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 
                 set.seed(58+q)
                 
                 ### required files:
                 source("~/functions.R")
                 source("skew normal copula.R")
                 source("C:/Users/u0125240/Documents/PhD/code/symmetry/symmetry tests.R")
                 source("C:/Users/u0125240/Documents/PhD/code/lineaire combinatie/Code part 1 Functions.R")
                 source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitLCS.R")
                 source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitPcopula.R")
                 source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitSPcopula.R")
                 source("C:/Users/u0125240/Documents/PhD/code/symmetry/fitNP.R")
                 

                 # Generate sample of pseudo-observations under the model with assumed central symmetry
                 MCU1.P <- rSNcopula(n = n1,alpha = rep(0,d1), rho = fit1P.S$par)
                 MCX1.P <- matrix(NA,nrow=n1,ncol=2)
                 MCX1.P[,1] <- qATD(beta = MCU1.P[,1],mu = margins1S.1$parameters[2]
                                    ,alpha = margins1S.1$parameters[1]
                                    ,phi = margins1S.1$parameters[3]
                                    ,nu = margins1S.1$parameters[4])
                 MCX1.P[,2] <- qATD(beta = MCU1.P[,2],mu = margins1S.2$parameters[2]
                                    ,alpha = margins1S.2$parameters[1]
                                    ,phi = margins1S.2$parameters[3]
                                    ,nu = margins1S.2$parameters[4])
                 
                 ### parametric fitting
                 # determining margins
                 MCmargins1R.1=marginfit(data = as.matrix(MCX1.P[,1]),crit = "AIC",symmetric=F)
                 MCmargins1R.2=marginfit(data = as.matrix(MCX1.P[,2]),crit = "AIC",symmetric=F)
                 
                 MCmargins1S.1=marginfit(data = as.matrix(MCX1.P[,1]),crit = "AIC",symmetric=T)
                 MCmargins1S.2=marginfit(data = as.matrix(MCX1.P[,2]),crit = "AIC",symmetric=T)
                 
                 # pseudo-observations
                 MCUP1.R = cbind(MCmargins1R.1$U,MCmargins1R.2$U)
                 MCUP1.S = cbind(MCmargins1S.1$U,MCmargins1S.2$U)
                 
                 # fits
                 MCfit1P.R = fitSNcopula(u = MCUP1.R,symmetric = F,nstart = 5,random = T) 
                 MCfit1P.S = fitSNcopula(u = MCUP1.S,symmetric = T,nstart = 5,random = T) 
                 
                 
                 ### semi-parametric fitting                 
                 # symmetrized sample w.r.t. the median
                 MCX1.NP <- jitter(X1AUG[sample(x = 1:(2*n1),size = n1,replace = T),],factor = 0.0001)
                 MCtheta1.NP <- apply(MCX1.NP,2,median)
                 MCX1.S <- sweep(x = -MCX1.NP,MARGIN = 2,STATS = 2*MCtheta1.NP,FUN = "+")
                 
                 MCU1.SP <- pobs(MCX1.NP)
                 MCU1.SPS = pobs(rbind(MCX1,MCX1.S))[1:n1,]
                 
                 
                 # regular
                 MCfit1SP.R = fitSNcopula(u = MCU1.SP,symmetric = F,nstart = 10,random = T)
                 # symmetric on symmetrized pseudo-observations
                 MCfit1SPS.S = fitSNcopula(u = MCU1.SPS,symmetric = T,nstart = 10,random = T)
                 
                 ### non-parametric fitting
                 MCfit1NP.R = fitNP(X = MCX1.NP,symmetric = F,IF = 2,GP = 100,ub = c(550,75),lb = c(50,15))
                 MCfit1NP.S = fitNP(X = MCX1.NP,symmetric = T,IF = 2,GP = 100,ub = c(550,75),lb = c(-200,0))
                 
                 
                 ### testing for symmetry
                 MCW1.P <- -2*(MCfit1P.S$ll-MCfit1P.R$ll)
                 MCW1.SPS <- -2*(MCfit1SPS.S$ll-MCfit1SP.R$ll)
                 MCW1.NP <- -2*(MCfit1NP.S$ll-MCfit1NP.R$ll)
                 MCW1.DHS <- DHS(X = MCX1.NP,type = "all",a = 2,b = 1,c = 1)
                 MCW1.DR <- DBST(X = MCX1.NP,depth = "H")
                 
                 out <- list("MCW1.P"=MCW1.P,"MCW1.SPS"=MCW1.SPS,"MCW1.NP"=MCW1.NP,"MCW1.DHS"=MCW1.DHS,"MCW1.DR"=MCW1.DR)
                 save(out,file=paste0("~/dataexamples/AsymptDistBodyMeasRun",q,".Rdata"))
               }
stopCluster(cl)



### Example 2: Log stock returns

library(quantmod)
getSymbols(c('^GDAXI','^FCHI','^N225','AMZN','NFLX','AMD'),from = "2010-01-01", to = "2018-12-31")
CAC40=as.numeric(weeklyReturn(FCHI,type="log"))[-470]
DAX=as.numeric(weeklyReturn(GDAXI,type="log"))
pairs(cbind(CAC40,DAX))
X2=as.matrix(cbind(DAX,CAC40))


### fitting a skew-normal copula with QBA-margins ###
#####################################################
seed=3158
set.seed(seed)

n2=nrow(X2)
d2=ncol(X2)

### parametric fitting
# determining margins
margins2R.1=marginfit(data = as.matrix(X2[,1]),crit = "AIC",symmetric=F) # t
margins2R.2=marginfit(data = as.matrix(X2[,2]),crit = "AIC",symmetric=F) # logistic

margins2S.1=marginfit(data = as.matrix(X2[,1]),crit = "AIC",symmetric=T) # t
margins2S.2=marginfit(data = as.matrix(X2[,2]),crit = "AIC",symmetric=T) # t

# pseudo-observations
UP2.R = cbind(margins2R.1$U,margins2R.2$U)
UP2.S = cbind(margins2S.1$U,margins2S.2$U)

# fits
fit2P.R = fitSNcopula(u = UP2.R,symmetric = F,nstart = 10,random = T) 
fit2P.S = fitSNcopula(u = UP2.S,symmetric = T,nstart = 10,random = T) 


### semi-parametric fitting

# symmetrized sample w.r.t. the median
theta2 <- apply(X2,2,median)
X2S <- sweep(x = -X2,MARGIN = 2,STATS = 2*theta2,FUN = "+")
X2AUG <- rbind(X2,X2S)

# pseudo-observations
USP2 = pobs(X2)
USPS2 = pobs(X2AUG)[1:n2,]

# regular
fit2SP.R = fitSNcopula(u = USP2,symmetric = F,nstart = 10,random = T)
# symmetric on regular pseudo-observations
fit2SP.S = fitSNcopula(u = USP2,symmetric = T,nstart = 10,random = T)
# symmetric on symmetrized pseudo-observations
fit2SPS.S = fitSNcopula(u = USPS2,symmetric = T,nstart = 10,random = T)

### non-parametric fitting
fit2NP.R = fitNP(X = X2,symmetric = F,IF = 2,GP = 100,ub = c(0.16,0.135),lb = c(-0.15,-0.13))
fit2NP.S = fitNP(X = X2,symmetric = T,IF = 2,GP = 100,ub = c(0.16,0.135),lb = c(-0.15,-0.13))


### testing for symmetry
W2.P <- -2*(fit2P.S$ll-fit2P.R$ll)
W2.SP <- -2*(fit2SP.S$ll-fit2SP.R$ll) # direct with chi-square
W2.SPS <- -2*(fit2SPS.S$ll-fit2SP.R$ll)
W2.NP <- -2*(fit2NP.S$ll-fit2NP.R$ll)
W2.DHS <- DHS(X = X2,type = "all",a = 2,b = 1,c = 1)
W2.DR <- DBST(X = X2+runif(n2*d2,-0.001,0.001),depth = "H")


### Asymptotic distribution of test statistics
dir.create("~/dataexamples")

testvalues=list("W2.P"=W2.P,"W2.SP"=W2.SP,"W2.SPS"=W2.SPS,"W2.NP"=W2.NP,"W2.DHS"=W2.DHS,"W2.DR"=W2.DR)
save(testvalues,file=paste0("~/dataexamples/StockTestValues.Rdata"))


# parametric
cl <- makeCluster(5)
registerDoParallel(cl)
### main loop for drawing sample from the model and refitting it to the sample
out <- foreach(q=1:400,.packages=c('sn','ks','copula','nloptr'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 
                 set.seed(58+q)
                 
                 ### required files:
                 source("~/functions.R")
                 
                 # Generate sample of pseudo-observations under the model with assumed central symmetry
                 MCU2.P <- rSNcopula(n = n2,alpha = rep(0,d2), rho = fit2P.S$par)
                 MCX2.P <- matrix(NA,nrow=n2,ncol=2)
                 MCX2.P[,1] <- qATD(beta = MCU2.P[,1],mu = margins2S.1$parameters[2]
                                    ,alpha = margins2S.1$parameters[1]
                                    ,phi = margins2S.1$parameters[3]
                                    ,nu = margins2S.1$parameters[4])
                 MCX2.P[,2] <- qATD(beta = MCU2.P[,2],mu = margins2S.2$parameters[2]
                                    ,alpha = margins2S.2$parameters[1]
                                    ,phi = margins2S.2$parameters[3]
                                    ,nu = margins2S.2$parameters[4])
                 
                 ### parametric fitting
                 # determining margins
                 MCmargins2R.1=marginfit(data = as.matrix(MCX2.P[,1]),crit = "AIC",symmetric=F)
                 MCmargins2R.2=marginfit(data = as.matrix(MCX2.P[,2]),crit = "AIC",symmetric=F)
                 
                 MCmargins2S.1=marginfit(data = as.matrix(MCX2.P[,1]),crit = "AIC",symmetric=T)
                 MCmargins2S.2=marginfit(data = as.matrix(MCX2.P[,2]),crit = "AIC",symmetric=T)
                 
                 # pseudo-observations
                 MCUP2.R = cbind(MCmargins2R.1$U,MCmargins2R.2$U)
                 MCUP2.S = cbind(MCmargins2S.1$U,MCmargins2S.2$U)
                 
                 # fits
                 MCfit2P.R = fitSNcopula(u = MCUP2.R,symmetric = F,nstart = 3,random = T) 
                 MCfit2P.S = fitSNcopula(u = MCUP2.S,symmetric = T,nstart = 3,random = T) 
                 
                 
                 ### semi-parametric fitting                 
                 # symmetrized sample w.r.t. the median
                 MCX2.NP <- jitter(X2AUG[sample(x = 1:(2*n2),size = n2,replace = T),],factor = 0.0001)
                 MCtheta2.NP <- apply(MCX2.NP,2,median)
                 MCX2.S <- sweep(x = -MCX2.NP,MARGIN = 2,STATS = 2*MCtheta2.NP,FUN = "+")
                 
                 MCU2.SP <- pobs(MCX2.NP)
                 MCU2.SPS = pobs(rbind(MCX2.NP,MCX2.S))[1:n2,]
                 
                 
                 # regular
                 MCfit2SP.R = fitSNcopula(u = MCU2.SP,symmetric = F,nstart = 10,random = T)
                 # symmetric on symmetrized pseudo-observations
                 MCfit2SPS.S = fitSNcopula(u = MCU2.SPS,symmetric = T,nstart = 10,random = T)
                 
                 ### non-parametric fitting
                 MCfit2NP.R = fitNP(X = MCX2.NP,symmetric = F,IF = 2,GP = 100,ub = c(550,75),lb = c(50,15))
                 MCfit2NP.S = fitNP(X = MCX2.NP,symmetric = T,IF = 2,GP = 100,ub = c(550,75),lb = c(-200,0))
                 
                 
                 ### testing for symmetry
                 MCW2.P <- -2*(MCfit2P.S$ll-MCfit2P.R$ll)
                 MCW2.SPS <- -2*(MCfit2SPS.S$ll-MCfit2SP.R$ll)
                 MCW2.NP <- -2*(MCfit2NP.S$ll-MCfit2NP.R$ll)
                 MCW2.DHS <- DHS(X = MCX2.NP,type = "all",a = 2,b = 1,c = 1)
                 MCW2.DR <- DBST(X = MCX2.NP,depth = "H")
                 
                 out <- list("MCW2.P"=MCW2.P,"MCW2.SPS"=MCW2.SPS,"MCW2.NP"=MCW2.NP,"MCW2.DHS"=MCW2.DHS,"MCW2.DR"=MCW2.DR)
                 save(out,file=paste0("~/dataexamples/AsymptDistStockRun",q,".Rdata"))
               }
stopCluster(cl)
