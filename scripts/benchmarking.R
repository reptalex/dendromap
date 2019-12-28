### denromap benchmarking
library(dendromap)
library(ggplot2)
library(mclust)


# Simulate dataset --------------------------------------------------------
set.seed(1)
m=500
n=20
row.tree <- rtree(m)
col.tree <- rtree(n)
S <- treeSim(5,row.tree,col.tree,row.depth.max = 2,
             row.depth.min=0.8,col.nodeset=c(21,23,24,33))

### count vs. presence-absence
ilogit <- function(x) 1/(1+exp(-x))
baseline_abundance=0
effect.size=10
noise_added=1
counts=1e4
compositional_link <- function(X,counts=1e4) apply(X,2,FUN=function(x,a) a*(exp(x)/sum(exp(x))),counts)
X <- (S$W %*% (effect.size*S$D) %*% t(S$V)) + rnorm(m*n,sd=noise_added)+baseline_abundance
N <- matrix(rpois(m*n,lambda=compositional_link(X,counts)),nrow=m,ncol=n)
# N <- matrix(rpois(m*n,lambda=exp(X)),nrow=m,ncol=n)
rownames(N) <- row.tree$tip.label
colnames(N) <- col.tree$tip.label

std <- function(x) (x-mean(x))/sd(x)
standardize <- function(X) apply(X,1,std) %>% apply(1,std)
N_std <- standardize(N)

# nb-like?  ---------------------------------------------------------------
mn <- rowMeans(N)
v <- apply(N,1,var)
plot(mn,v,log='xy')
abline(0,2)

## ish...

# U distn and Pvalues -----------------------------------------------------

rc <- makeRCtable(N,row.tree,col.tree)
rc[,rc_ix:=paste(row.node,col.node,sep='_')]
s <- S$Lineages
s[,rc_ix:=paste(row.node,col.node,sep='_')]

rc[,label:='null']
rc[rc_ix %in% s$rc_ix,label:='real']

ggplot(rc,aes(P,fill=label))+geom_density(alpha=0.6)
# 
# W=treeBasis(row.tree)
# V=treeBasis(col.tree)
# U <- t(W) %*% N %*% V
# U_null <- t(W) %*% N[sample(m),sample(n)] %*% V
# D <- data.table('u'=log(abs(c(c(U),c(U_null)))),
#                 'dataset'=rep(c('real','null'),each=ncol(W)*ncol(V)))
# 
# D[u<=-10,u:=NA]
# 
# u <- D[dataset=='real' & !is.na(u),u]
# md=median(u)
# ### can we estimate a Gaussian based on the ecdf below the median?
# 
# hist(y,breaks=100)
# abline(v=md,lwd=5,col=rgb(0,0.5,0.2))
# 
# ### zero-censored normal fitting for y?
# # y = rnorm(10000,0.5,2) # simulate some sample data 
# y=md-u
# y[y<0]=0 # censor at zero
# minlogL = function (p ) { 
#   mu = p[1]
#   logsigma=p[2] # must used log-sigma because the second parameter can be negative while sigma must be positive
#   return(sum(-log(dnorm(y[y>0],mu,exp(logsigma)))) - sum(y<=0)*log(pnorm(0,mu,exp(logsigma))))
# }
# fit=nlm(minlogL,c(0,0)) 
# mn=-fit$estimate[1]
# sig2=exp(fit$estimate[2])
# 
# xx=seq(min(u),max(u),length.out = 1e3)
# hist(u,freq = F,breaks=1e2)
# lines(xx,dnorm(xx,mn+md,sig2),lwd=2,col=rgb(0,0.5,0.2))
# 
# 
# null_cdf <-  ecdf(D[dataset=='null']$u)
# D[,P:=1-null_cdf(u)]
# 
# 
# 
# if (any(D$P==0)){
#   nn=sum(D$P==0)
#   base::cat(paste('\n',nn,' P-values were 0. Will estimate tail probabilities assuming log(u^2)~rnorm(mu,sd). You may want to consider manually increasing n_sim for more accurate null distribution',sep=''))
#   y <- y[y>-20]
#   mu <- mean(y)
#   sig <- sd(y)
#   min.P <- D[P>0,min(P)]
#   u.min <- max(D[P>0][P==min(P),u])
#   estimate_tail <- function(y1,mu,sig,ymin,pmin) (1-pnorm(y1,mu,sig))/(1-pnorm(ymin,mu,sig))*pmin
#   D[P==0,P:=estimate_tail(u,mu,sig,u.min,min.P)]
# }
# 
# 
# 
# ggplot(D,aes(u,fill=dataset))+
#   geom_histogram(alpha=0.3,position = 'identity',bins = 100)
# 
# f <- ecdf(c(U_null))


# dendromap ---------------------------------------------------------------


x <- dendromap(N,row.tree,col.tree,Pval_threshold = 0.1)


y <- S$W %*% S$D %*% t(S$V)
y[y==0] <- NA
amplifier=3
plot.dendromap(S,ilogit(amplifier*y),orient.nodes=F)
dev.new()
plot.dendromap(x,orient.nodes=F)
