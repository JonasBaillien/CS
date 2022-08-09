library(ggplot2)
library(gridExtra)

### script of plotting bivariate skew-normal copula ###
#######################################################

n=2000
u=data.frame(rSNcopula(n = n,alpha = c(0,0),rho = -0.8))
colnames(u)=c("x","y")
p1 <- ggplot(data = u, aes(x = x,y = y)) + 
  geom_point() +
  ylab(expression(u[2])) +
  xlab(expression(u[1])) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14))

u2=data.frame(rSNcopula(n = n,alpha = c(0,0),rho = 0.0))
colnames(u2)=c("x","y")
p2 <- ggplot(data = u2, aes(x = x,y = y)) + 
  geom_point() +
  ylab(expression(u[2])) +
  xlab(expression(u[1])) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14))

u3=data.frame(rSNcopula(n = n,alpha = c(0,0),rho = 0.8))
colnames(u3)=c("x","y")
p3 <- ggplot(data = u3, aes(x = x,y = y)) + 
  geom_point() +
  ylab(expression(u[2])) +
  xlab(expression(u[1])) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14))

u4=data.frame(rSNcopula(n = n,alpha = c(5,5),rho = -0.8))
colnames(u4)=c("x","y")
p4 <- ggplot(data = u4, aes(x = x,y = y)) + 
  geom_point() +
  ylab(expression(u[2])) +
  xlab(expression(u[1])) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14))

u5=data.frame(rSNcopula(n = n,alpha = c(5,5),rho = 0.0))
colnames(u5)=c("x","y")
p5 <- ggplot(data = u5, aes(x = x,y = y)) + 
  geom_point() +
  ylab(expression(u[2])) +
  xlab(expression(u[1])) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14))

u6=data.frame(rSNcopula(n = n,alpha = c(5,5),rho = 0.8))
colnames(u6)=c("x","y")
p6 <- ggplot(data = u6, aes(x = x,y = y)) + 
  geom_point() +
  ylab(expression(u[2])) +
  xlab(expression(u[1])) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14))

u7=data.frame(rSNcopula(n = n,alpha = c(-5,-5),rho = -0.8))
colnames(u7)=c("x","y")
p7 <- ggplot(data = u7, aes(x = x,y = y)) + 
  geom_point() +
  ylab(expression(u[2])) +
  xlab(expression(u[1])) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14))

u8=data.frame(rSNcopula(n = n,alpha = c(-5,-5),rho = 0.0))
colnames(u8)=c("x","y")
p8 <- ggplot(data = u8, aes(x = x,y = y)) + 
  geom_point() +
  ylab(expression(u[2])) +
  xlab(expression(u[1])) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14))

u9=data.frame(rSNcopula(n = n,alpha = c(-5,-5),rho = 0.8))
colnames(u9)=c("x","y")
p9 <- ggplot(data = u9, aes(x = x,y = y)) + 
  geom_point() +
  ylab(expression(u[2])) +
  xlab(expression(u[1])) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14))


u10=data.frame(rSNcopula(n = n,alpha = c(-5,5),rho = -0.8))
colnames(u10)=c("x","y")
p10 <- ggplot(data = u10, aes(x = x,y = y)) + 
  geom_point() +
  ylab(expression(u[2])) +
  xlab(expression(u[1])) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14))

u11=data.frame(rSNcopula(n = n,alpha = c(-5,5),rho = 0.0))
colnames(u11)=c("x","y")
p11 <- ggplot(data = u11, aes(x = x,y = y)) + 
  geom_point() +
  ylab(expression(u[2])) +
  xlab(expression(u[1])) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14))

u12=data.frame(rSNcopula(n = n,alpha = c(-5,5),rho = 0.8))
colnames(u12)=c("x","y")
p12 <- ggplot(data = u12, aes(x = x,y = y)) + 
  geom_point() +
  ylab(expression(u[2])) +
  xlab(expression(u[1])) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14))

u13=data.frame(rSNcopula(n = n,alpha = c(5,-5),rho = -0.8))
colnames(u13)=c("x","y")
p13 <- ggplot(data = u13, aes(x = x,y = y)) + 
  geom_point() +
  ylab(expression(u[2])) +
  xlab(expression(u[1])) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14))

u14=data.frame(rSNcopula(n = n,alpha = c(5,-5),rho = 0.0))
colnames(u14)=c("x","y")
p14 <- ggplot(data = u14, aes(x = x,y = y)) + 
  geom_point() +
  ylab(expression(u[2])) +
  xlab(expression(u[1])) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14))

u15=data.frame(rSNcopula(n = n,alpha = c(5,-5),rho = 0.8))
colnames(u15)=c("x","y")
p15 <- ggplot(data = u15, aes(x = x,y = y)) + 
  geom_point() +
  ylab(expression(u[2])) +
  xlab(expression(u[1])) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14))

x11()
grid.arrange(p7,p10,p1,p13,p4,p8,p11,p2,p14,p5,p9,p12,p3,p15,p6, nrow = 3,ncol = 5)
