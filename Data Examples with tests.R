### load in required functions ###
##################################
source("functions.R")


#### Data examples ###
######################

### Example 1: Body measeruments data
# loading data (from lessR package)
library(lessR)
data("dataBodyMeas")
body <- dataBodyMeas

# plotting data
x11()
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

### non-parametric fitting
# unconstrained
fit1NP.R = fitNP(X = X1,symmetric = F,IF = 2,GP = 100)
# centrally symmetric
fit1NP.S = fitNP(X = X1,symmetric = T,IF = 2,GP = 100)

## symmetrized sample w.r.t. the mode
theta1 <- fit1NP.S$mode
X1S <- sweep(x = -X1,MARGIN = 2,STATS = 2*theta1,FUN = "+")
X1AUG <- rbind(X1,X1S)

### semi-parametric fitting
## pseudo-observations
# empirical original data
USP1 = pobs(X1)
# empirical in augmented data
USPS1 = pobs(X1AUG)[1:n1,]

## fitting SN-copula
# regular
fit1SP1.R = fitSNcopula(u = USP1,symmetric = F,nstart = 10,random = T)
# symmetric on symmetrized pseudo-observations
fit1SP1.S = fitSNcopula(u = USPS1,symmetric = T,nstart = 10,random = T)

## fitting ST-copula
# regular
fit1SP2.R = fitSTcopula(u = USP1,symmetric = F,nstart = 10,random = T)
# symmetric on regular pseudo-observations
fit1SP2.S = fitSTcopula(u = USPS1,symmetric = T,nstart = 10,random = T)

### testing for symmetry
W1.P <- -2*(fit1P.S$ll-fit1P.R$ll)
W1.SP1 <- -2*(fit1SP1.S$ll-fit1SP1.R$ll) # direct with chi-square
W1.SP2 <- -2*(fit1SP2.S$ll-fit1SP2.R$ll) # direct with chi-square
W1.NP <- -2*(fit1NP.S$ll-fit1NP.R$ll)
W1.DHS <- DHS(X = X1,type = "all",a = 2,b = 1,c = 1)
W1.DR <- DBST(X = X1+runif(n1*d1,-0.001,0.001),depth = "H")

### Asymptotic distribution of test statistics
dir.create("dataexamples")

testvalues=list("W1.P"=W1.P,"W1.SP1"=W1.SP1,"W1.SP2"=W1.SP2,"W1.NP"=W1.NP,"W1.DHS"=W1.DHS,"W1.DR"=W1.DR)
save(testvalues,file=paste0("dataexamples/bodymeas/BodyMeasTestValues.Rdata"))


# parametric tests in parallel
cl <- makeCluster(5)
registerDoParallel(cl)
### main loop for drawing sample from the model and refitting it to the sample
out <- foreach(q=1:400,.packages=c('sn','ks','copula','nloptr'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 
                 set.seed(58+q)
                 
                 ### required files:
                 source("functions.R")
                 
                 
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
                 
                 
                 
                 
                 ### non-parametric fitting
                 MCX1.NP <- jitter(X1AUG[sample(x = 1:(2*n1),size = n1,replace = T),],factor = 0.0001)
                 MCfit1NP.R = fitNP(X = MCX1.NP,symmetric = F,IF = 2,GP = 100)
                 MCfit1NP.S = fitNP(X = MCX1.NP,symmetric = T,IF = 2,GP = 100)
                 
                 
                 # symmetrized sample w.r.t. the mode
                 MCtheta1.NP <- MCfit1NP.S$mode
                 MCX1.S <- sweep(x = -MCX1.NP,MARGIN = 2,STATS = 2*MCtheta1.NP,FUN = "+")
                 
                 MCU1.SP <- pobs(MCX1.NP)
                 MCU1.SPS = pobs(rbind(MCX1,MCX1.S))[1:n1,]
                 
                 ### semi-parametric fitting 
                 # regular
                 MCfit1SP1.R = fitSNcopula(u = MCU1.SP,symmetric = F,nstart = 10,random = T)
                 # symmetric on symmetrized pseudo-observations
                 MCfit1SP1.S = fitSNcopula(u = MCU1.SPS,symmetric = T,nstart = 10,random = T)
                 # regular
                 MCfit1SP2.R = fitSTcopula(u = MCU1.SP,symmetric = F,nstart = 10,random = T)
                 # symmetric on symmetrized pseudo-observations
                 MCfit1SP2.S = fitSTcopula(u = MCU1.SPS,symmetric = T,nstart = 10,random = T)
                 
                 ### testing for symmetry
                 MCW1.P <- -2*(MCfit1P.S$ll-MCfit1P.R$ll)
                 MCW1.SP1 <- -2*(MCfit1SP1.S$ll-MCfit1SP1.R$ll)
                 MCW1.SP2 <- -2*(MCfit1SP2.S$ll-MCfit1SP2.R$ll)
                 MCW1.NP <- -2*(MCfit1NP.S$ll-MCfit1NP.R$ll)
                 MCW1.DHS <- DHS(X = MCX1.NP,type = "all",a = 2,b = 1,c = 1)
                 MCW1.DR <- DBST(X = MCX1.NP,depth = "H")
                 
                 out <- list("MCW1.P"=MCW1.P,"MCW1.SP1"=MCW1.SP1,"MCW1.SP2"=MCW1.SP2,"MCW1.NP"=MCW1.NP,"MCW1.DHS"=MCW1.DHS,"MCW1.DR"=MCW1.DR)
                 save(out,file=paste0("dataexamples/AsymptDistBodyMeasRun",q,".Rdata"))
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

### non-parametric fitting
fit2NP.R = fitNP2(X = X2,symmetric = F,IF = 2,GP = 100)
fit2NP.S = fitNP2(X = X2,symmetric = T,IF = 2,GP = 100)
# symmetrized sample w.r.t. the mode
theta2 <- fit2NP.S$mode
X2S <- sweep(x = -X2,MARGIN = 2,STATS = 2*theta2,FUN = "+")
X2AUG <- rbind(X2,X2S)
# symmetrized sample w.r.t. the median
theta2 <- apply(X2,2,median)
X2S <- sweep(x = -X2,MARGIN = 2,STATS = 2*theta2,FUN = "+")
X2AUG <- rbind(X2,X2S)

### semi-parametric fitting
# pseudo-observations
USP2 = pobs(X2)
USPS2 = pobs(X2AUG)[1:n2,]

# regular
fit2SP1.R = fitSNcopula(u = USP2,symmetric = F,nstart = 10,random = T)
# symmetric on regular pseudo-observations
fit2SP1.S = fitSNcopula(u = USPS2,symmetric = T,nstart = 10,random = T)

# regular
fit2SP2.R = fitSTcopula(u = USP2,symmetric = F,nstart = 10,random = T)
# symmetric on regular pseudo-observations
fit2SP2.S = fitSTcopula(u = USPS2,symmetric = T,nstart = 10,random = T)

### testing for symmetry
W2.P <- -2*(fit2P.S$ll-fit2P.R$ll)
W2.SP1 <- -2*(fit2SP1.S$ll-fit2SP1.R$ll) 
W2.SP2 <- -2*(fit2SP2.S$ll-fit2SP2.R$ll) 
W2.NP <- -2*(fit2NP.S$ll-fit2NP.R$ll)
W2.DHS <- DHS(X = X2,type = "all",a = 2,b = 1,c = 1)
W2.DR <- DBST(X = X2+runif(n2*d2,-0.001,0.001),depth = "H")


### Asymptotic distribution of test statistics
dir.create("dataexamples")

testvalues=list("W2.P"=W2.P,"W2.SP1"=W2.SP1,"W2.SP2"=W2.SP2,"W2.NP"=W2.NP,"W2.DHS"=W2.DHS,"W2.DR"=W2.DR)
save(testvalues,file="~/PhD/code/symmetry/output/stock/StockTestValues.Rdata")


# parametric
cl <- makeCluster(5)
registerDoParallel(cl)
### main loop for drawing sample from the model and refitting it to the sample
out <- foreach(q=1:400,.packages=c('sn','ks','copula','nloptr'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 
                 set.seed(58+q)
                 
                 ### required files:
                 source("functions.R")
                 
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
                 
                 
                 
                 ### non-parametric fitting
                 MCX2.NP <- jitter(X2AUG[sample(x = 1:(2*n2),size = n2,replace = T),],factor = 0.0001)
                 MCfit2NP.R = fitNP(X = MCX2.NP,symmetric = F,IF = 2,GP = 100)
                 MCfit2NP.S = fitNP(X = MCX2.NP,symmetric = T,IF = 2,GP = 100)
                 
                 # symmetrized sample w.r.t. the mode
                 MCtheta2.NP <- MCfit2NP.S$mode
                 MCX2.S <- sweep(x = -MCX2.NP,MARGIN = 2,STATS = 2*MCtheta2.NP,FUN = "+")
                 
                 MCU2.SP <- pobs(MCX2.NP)
                 MCU2.SPS = pobs(rbind(MCX2.NP,MCX2.S))[1:n2,]
                 
                 ### semi-parametric fitting  
                 # regular
                 MCfit2SP1.R = fitSNcopula(u = MCU2.SP,symmetric = F,nstart = 10,random = T)
                 # symmetric on symmetrized pseudo-observations
                 MCfit2SP1.S = fitSNcopula(u = MCU2.SPS,symmetric = T,nstart = 10,random = T)
                 # regular
                 MCfit2SP2.R = fitSTcopula(u = MCU2.SP,symmetric = F,nstart = 10,random = T)
                 # symmetric on symmetrized pseudo-observations
                 MCfit2SP2.S = fitSTcopula(u = MCU2.SPS,symmetric = T,nstart = 10,random = T)
                 
                 ### testing for symmetry
                 MCW2.P <- -2*(MCfit2P.S$ll-MCfit2P.R$ll)
                 MCW2.SP1 <- -2*(MCfit2SP1.S$ll-MCfit2SP1.R$ll)
                 MCW2.SP2 <- -2*(MCfit2SP2.S$ll-MCfit2SP2.R$ll)
                 MCW2.NP <- -2*(MCfit2NP.S$ll-MCfit2NP.R$ll)
                 MCW2.DHS <- DHS(X = MCX2.NP,type = "all",a = 2,b = 1,c = 1)
                 MCW2.DR <- DBST(X = MCX2.NP,depth = "H")
                 
                 out <- list("MCW2.P"=MCW2.P,"MCW2.SP1"=MCW2.SP1,"MCW2.SP2"=MCW2.SP2,"MCW2.NP"=MCW2.NP,"MCW2.DHS"=MCW2.DHS,"MCW2.DR"=MCW2.DR)
                 save(out,file=paste0("dataexamples/AsymptDistStockRun",q,".Rdata"))
               }
stopCluster(cl)











### Example 3: Wine data
# loading data (from lessR package)
winequality.red <- read.csv("~/PhD/code/datasets/winequality-red.csv", sep=";")
winequality.white <- read.csv("~/PhD/code/datasets/winequality-white.csv", sep=";")

# white wine of quality 4
ww4=winequality.white[winequality.white$quality==4,c(2,9,10)] 
# red wine of quality 4
rw4=winequality.red[winequality.red$quality==4,c(2,9,10)] 
# white wine of quality 7
ww7=winequality.white[winequality.white$quality==7,c(2,9,10)] 
# red wine of quality 7
rw7=winequality.red[winequality.red$quality==7,c(2,9,10)] 

# data
X3=as.matrix(ww7)   

# dimensions of data
d3=length(X3[1,])
n3=length(X3[,1])

### load in required functions ###
##################################
source("functions.R")



### fitting a skew-normal copula with QBA-margins ###
#####################################################
seed=3158
set.seed(seed)

n3=nrow(X3)
d3=ncol(X3)

### parametric fitting

## determining marginals
# unconstrained
margins3R.1=marginfit(data = as.matrix(X3[,1]),crit = "AIC",symmetric=F) # normal
margins3R.2=marginfit(data = as.matrix(X3[,2]),crit = "AIC",symmetric=F) # normal
margins3R.3=marginfit(data = as.matrix(X3[,3]),crit = "AIC",symmetric=F) # logistic


# symmetric
margins3S.1=marginfit(data = as.matrix(X3[,1]),crit = "AIC",symmetric=T) # t
margins3S.2=marginfit(data = as.matrix(X3[,2]),crit = "AIC",symmetric=T) # normal
margins3S.3=marginfit(data = as.matrix(X3[,3]),crit = "AIC",symmetric=T) # logistic

## pseudo-observations
# unconstrained
UP3.R = cbind(margins3R.1$U,margins3R.2$U,margins3R.3$U)
# symmetric
UP3.S = cbind(margins3S.1$U,margins3S.2$U,margins3S.3$U)


## fits
# unconstrained
fit3P.R = fitSNcopula(u = UP3.R,symmetric = F,nstart = 10,random = T) 
# symmetric
fit3P.S = fitSNcopula(u = UP3.S,symmetric = T,nstart = 10,random = T) 


### non-parametric fitting
# unconstrained
fit3NP.R = fitNP(X = X3,symmetric = F,IF = 2,GP = 100)
# centrally symmetric
fit3NP.S = fitNP(X = X3,symmetric = T,IF = 2,GP = 100)

theta3 <- fit3NP.S$mode
X3S <- sweep(x = -X3,MARGIN = 2,STATS = 2*theta3,FUN = "+")
X3AUG <- rbind(X3,X3S)

## pseudo-observations
# empirical original data
USP3 = pobs(X3)
# empirical in augmented data
USPS3 = pobs(X3AUG)[1:n3,]

## fitting SN-copula
# regular
fit3SP1.R = fitSNcopula(u = USP3,symmetric = F,nstart = 5,random = T)
# symmetric on symmetrized pseudo-observations
fit3SP1.S = fitSNcopula(u = USPS3,symmetric = T,nstart = 20,random = T)

## fitting ST-copula
# regular
fit3SP2.R = fitSTcopula(u = USP3,symmetric = F,nstart = 5,random = T)
# symmetric on regular pseudo-observations
fit3SP2.S = fitSTcopula(u = USPS3,symmetric = T,nstart = 20,random = T)


### testing for symmetry
W3.P <- -2*(fit3P.S$ll-fit3P.R$ll)
W3.SP1 <- -2*(fit3SP1.S$ll-fit3SP1.R$ll) # direct with chi-square
W3.SP2 <- -2*(fit3SP2.S$ll-fit3SP2.R$ll) # direct with chi-square
W3.NP <- -2*(fit3NP.S$ll-fit3NP.R$ll)
W3.DHS <- DHS(X = as.matrix(X3),type = "all",a = 2,b = 1,c = 1)
W3.DR <- DBST(X = as.matrix(X3)+runif(n3*d3,-0.001,0.001),depth = "H")

### Asymptotic distribution of test statistics
dir.create("dataexamples")

testvalues=list("W3.P"=W3.P,"W3.SP1"=W3.SP1,"W3.SP2"=W3.SP2,"W1.NP"=W1.NP,"W1.DHS"=W1.DHS,"W1.DR"=W1.DR)
save(testvalues,file=paste0("dataexamples/WW7.Rdata"))

# parametric tests in parallel
cl <- makeCluster(5)
registerDoParallel(cl)
### main loop for drawing sample from the model and refitting it to the sample
out <- foreach(q=indmis,.packages=c('sn','ks','copula','nloptr'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 
                 set.seed(58000+q)
                 
                 ### required files:
                 source("functions.R")
                 
                 
                 
                 # Generate sample of pseudo-observations under the model with assumed central symmetry
                 MCU3.P <- rSNcopula(n = n3,alpha = rep(0,d3), rho = fit3P.S$par)
                 MCX3.P <- matrix(NA,nrow=n3,ncol=d3)
                 MCX3.P[,1] <- qATD(beta = MCU3.P[,1],mu = margins3S.1$parameters[2]
                                    ,alpha = margins3S.1$parameters[1]
                                    ,phi = margins3S.1$parameters[3]
                                    ,nu = margins3S.1$parameters[4])
                 MCX3.P[,2] <- qAND(beta = MCU3.P[,2],mu = margins3S.2$parameters[2]
                                    ,alpha = margins3S.2$parameters[1]
                                    ,phi = margins3S.2$parameters[3])
                 MCX3.P[,3] <- qALoD(beta = MCU3.P[,3],mu = margins3S.3$parameters[2]
                                     ,alpha = margins3S.3$parameters[1]
                                     ,phi = margins3S.3$parameters[3])
                 
                 ### parametric fitting
                 # determining margins
                 MCmargins3R.1=marginfit(data = as.matrix(MCX3.P[,1]),crit = "AIC",symmetric=F)
                 MCmargins3R.2=marginfit(data = as.matrix(MCX3.P[,2]),crit = "AIC",symmetric=F)
                 MCmargins3R.3=marginfit(data = as.matrix(MCX3.P[,3]),crit = "AIC",symmetric=F)
                 
                 MCmargins3S.1=marginfit(data = as.matrix(MCX3.P[,1]),crit = "AIC",symmetric=T)
                 MCmargins3S.2=marginfit(data = as.matrix(MCX3.P[,2]),crit = "AIC",symmetric=T)
                 MCmargins3S.3=marginfit(data = as.matrix(MCX3.P[,3]),crit = "AIC",symmetric=T)
                 
                 
                 # pseudo-observations
                 MCUP3.R = cbind(MCmargins3R.1$U,MCmargins3R.2$U,MCmargins3R.3$U)
                 MCUP3.S = cbind(MCmargins3S.1$U,MCmargins3S.2$U,MCmargins3S.3$U)
                 
                 # fits
                 MCfit3P.R = fitSNcopula(u = MCUP3.R,symmetric = F,nstart = 5,random = T) 
                 MCfit3P.S = fitSNcopula(u = MCUP3.S,symmetric = T,nstart = 5,random = T) 
                 
                 ### non-parametric fitting
                 MCX3.NP <- jitter(X3AUG[sample(x = 1:(2*n3),size = n3,replace = T),],factor = 0.0001)
                 MCfit3NP.R = fitNP(X = MCX3.NP,symmetric = F,IF = 2,GP = 100)
                 MCfit3NP.S = fitNP(X = MCX3.NP,symmetric = T,IF = 2,GP = 100)
                 
                 # symmetrized sample w.r.t. the mode
                 MCtheta3.NP <- MCfit3NP.S$mode
                 MCX3.S <- sweep(x = -MCX3.NP,MARGIN = 2,STATS = 2*MCtheta3.NP,FUN = "+")
                 
                 
                 ### semi-parametric fitting 
                 MCU3.SP <- pobs(MCX3.NP)
                 MCU3.SPS = pobs(rbind(MCX3.NP,MCX3.S))[1:n3,]
        
                 # regular fitting
                 fitReg.SP1 <- fitSNcopula(u = MCU3.SP,symmetric = F,nstart = 7,random=T)
                 # symmetric fitting with assumed symmetric margins
                 fitSym.SP1 <- fitSNcopula(u = MCU3.SPS,symmetric = T,nstart = 20,random=T)
                 # regular fitting
                 fitReg.SP2 <- fitSTcopula(u = MCU3.SP,symmetric = F,nstart = 2,random = T)
                 # symmetric fitting with assumed symmetric margins
                 fitSym.SP2 <- fitSTcopula(u = MCU3.SPS,symmetric = T,nstart = 20,random = T)
                 
                 
                 
                 ### testing for symmetry
                 MCW3.P <- -2*(MCfit3P.S$ll-MCfit3P.R$ll)
                 MCW3.SP1 <- -2*(fitSym.SP1$ll-fitReg.SP1$ll)  
                 MCW3.SP2 <- -2*(fitSym.SP2$ll-fitReg.SP2$ll)  
                 MCW3.NP <- -2*(MCfit3NP.S$ll-MCfit3NP.R$ll)
                 MCW3.DHS <- DHS(X = MCX3.NP,type = "all",a = 2,b = 1,c = 1)
                 MCW3.DR <- DBST(X = MCX3.NP,depth = "H")
                 
                 out <- list("MCW3.P"=MCW3.P,"MCW3.SP1"=MCW3.SP1,"MCW3.SP2"=MCW3.SP2,"MCW3.NP"=MCW3.NP,"MCW3.DHS"=MCW3.DHS,"MCW3.DR"=MCW3.DR)
                 save(out,file=paste0("dataexamples/WW7_",q,".Rdata"))
               }
stopCluster(cl)





