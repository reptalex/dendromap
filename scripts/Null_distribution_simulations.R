library(phylofactor)
library(ggplot2)
# library(dendromap) ## must be loaded atm

rData <- function(m=1e3,n=10,mean=0,sd=1,mean_counts=1e4,presence.absence=F){
  clr_inv <- function(x) exp(x)/sum(exp(x))
  X <- matrix(rnorm(n=m*n,mean=mean,sd=sd),nrow=m) %>%
    apply(2,clr_inv) %>% 
    apply(2,FUN=function(p) rmultinom(1,rpois(1,lambda=mean_counts),p))
  if (presence.absence){
    X[X>0] <- 1
  }
  return(X)
}



# Approximate log-normality -----------------------------------------------
set.seed(1)
m=1e3
n=1e2
row.tree <- rtree(m)
col.tree <- rtree(n)

V <- treeBasis(row.tree)
W <- treeBasis(col.tree)

mu=rnorm(m,sd=1)
X <- rData(m,n,mean=mu,sd=3,presence.absence=F)

mns <- rowMeans(X)
vs <- apply(X,1,var)
plot(mns,vs,xlab='mean',ylab='var',log='xy')
abline(0,1)
abline(1,2,col='blue',lwd=2)
glm(log(vs)~log(mns)) ## slope >2


plot(sort(rowSums(X),decreasing = T),log='y')


U <- t(V) %*% X %*% W
y <- c(U)
z <- log(y^2)[log(y^2)>-10]
x <- seq(min(z),max(z),length.out=1e3)
png(filename = 'Figures/Null_distribution_of_U.png',height = 7,width=10,units='in',res = 200)
  par(mfrow=c(1,2))
  hist(z,breaks=100,freq = F,xlab='log(u^2)')
  lines(x,dnorm(x,mean=mean(z),sd=sd(z)),col='blue')
  qqplot(z,rnorm(1e3,mean=mean(z),sd=sd(z)),xlab='log(u^2)',ylab='Gaussian Quantiles')
  abline(0,1)
dev.off()

var(c(U))
sqrt(m*n)*(0.1*(1-0.1))  ## sqrt(m*n)*var(X[i,j])


# presence-absence data ---------------------------------------------------
mu=rnorm(m,sd=0.2)
X <- rData(m,n,mean=mu,sd=1,presence.absence=T)

mns <- rowMeans(X)
vs <- apply(X,1,var)
plot(mns,vs,xlab='mean',ylab='var',log='xy')
abline(0,1)
abline(1,2,col='blue',lwd=2)
glm(log(vs)~log(mns)) ## slope >2


plot(sort(rowSums(X),decreasing = T),log='y')


U <- t(V) %*% X %*% W
y <- c(U)
y <- y[y!=0]
z <- log(y^2)[log(y^2)>-10]
x <- seq(min(z),max(z),length.out=1e3)
# png(filename = 'Figures/Null_distribution_of_U.png',height = 7,width=10,units='in',res = 200)
par(mfrow=c(1,2))
hist(z,breaks=30,freq = F,xlab='log(u^2)')
lines(x,dnorm(x,mean=mean(z),sd=sd(z)),col='blue')
qqplot(z,rnorm(1e3,mean=mean(z),sd=sd(z)),xlab='log(u^2)',ylab='Gaussian Quantiles')
abline(0,1)
# dev.off()

var(c(U))
sqrt(m*n)*(0.1*(1-0.1))  ## sqrt(m*n)*var(X[i,j])



# mean-variance relationships ---------------------------------------------


mset <- round(10^seq(1,3,length.out=7))
nset <- 2^(2:6)

DF <- expand.grid('m'=mset,'n'=nset)
DF$mean <- 0
DF$sd <- 0

for (i in 1:nrow(DF)){
  
  m <- DF$m[i]
  n <- DF$n[i]
  
  mu=rnorm(m,sd=1)
  X <- rData(m,n,mean=mu,sd=3,presence.absence=F)
  V <- treeBasis(rtree(m))
  W <- treeBasis(rtree(n))
  U <- t(V) %*% X %*% W
  y <- c(U)
  z <- log(y^2)[log(y^2)>-10]
  DF$mean[i] <- mean(z)
  DF$sd[i] <- sd(z)
}

names(DF)[2] <- "nContinents"
ggplot(DF,aes(m,mean,by=nContinents))+
  geom_point()+
  geom_smooth(method='glm')+
  facet_grid(.~nContinents)+
  scale_x_continuous('Number of Species',trans='log',breaks = 5^(0:4))


ggplot(DF,aes(m,mean,color=n))+
  geom_point()+
  geom_smooth(method='glm')+
  scale_x_continuous(trans='log',breaks = 5^(0:4))

## the mean shrinks with log(m)

glm(mean~log(m),data=DF)
## slope = -2.2

ggplot(DF,aes(m,sd))+
  geom_point()+
  geom_smooth(method='glm')+
  scale_x_continuous(trans='log',breaks = 5^(0:4))



# Fingerprint of true positives -------------------------------------------
rm(list=ls())
gc()
rData <- function(m=1e3,n=10,mean=0,sd=1,mean_counts=1e4,presence.absence=F){
  clr_inv <- function(x) exp(x)/sum(exp(x))
  X <- matrix(rnorm(n=m*n,mean=mean,sd=sd),nrow=m) %>%
    apply(2,clr_inv) %>% 
    apply(2,FUN=function(p) rmultinom(1,rpois(1,lambda=mean_counts),p))
  if (presence.absence){
    X[X>0] <- 1
  }
  return(X)
}


m <- 1e3
n <- 20
set.seed(1)
row.tree <- rtree(m)
col.tree <- rtree(n)
########################### SIM 
S <- treeSim(20,row.tree,col.tree)




############
U <- t(S$W) %*% rData(m,n,mean_counts=5e3) %*% S$V
y <- c(U)
y <- y[y!=0]
z <- log(y^2)[log(y^2)>-10]

DF <- data.frame('z'=z,'effect_size'=0)

X <- S$W %*% (4*S$D) %*% t(S$V)
X <- X+rnorm(m*n,sd=1)
clr_inv <- function(x) exp(x)/sum(exp(x))
X <- apply(X,2,clr_inv) %>% 
        apply(2,FUN=function(p) rmultinom(1,rpois(1,lambda=5e3),p))

U <- t(S$W) %*% X %*% S$V
y <- c(U)
y <- y[y!=0]
z <- log(y^2)[log(y^2)>-10]
DF <- rbind(DF,data.frame('z'=z,'effect_size'=4))


X <- S$W %*% (8*S$D) %*% t(S$V)
X <- X+rnorm(m*n,sd=1)
clr_inv <- function(x) exp(x)/sum(exp(x))
X <- apply(X,2,clr_inv) %>% 
  apply(2,FUN=function(p) rmultinom(1,rpois(1,lambda=5e3),p))

U <- t(S$W) %*% X %*% S$V
y <- c(U)
y <- y[y!=0]
z <- log(y^2)[log(y^2)>-10]
DF <- rbind(DF,data.frame('z'=z,'effect_size'=8))


X <- S$W %*% (16*S$D) %*% t(S$V)
X <- X+rnorm(m*n,sd=1)
clr_inv <- function(x) exp(x)/sum(exp(x))
X <- apply(X,2,clr_inv) %>% 
  apply(2,FUN=function(p) rmultinom(1,rpois(1,lambda=5e3),p))

U <- t(S$W) %*% X %*% S$V
y <- c(U)
y <- y[y!=0]
z <- log(y^2)[log(y^2)>-10]
DF <- rbind(DF,data.frame('z'=z,'effect_size'=16))
DF <- as.data.table(DF)

ggplot(DF,aes(z,group=effect_size,color=effect_size,fill=effect_size))+
  geom_density(alpha=0.5)+
  ggtitle('Impact of True-Positives on Distribution of test-statistic')




