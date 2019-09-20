library(dendromap)
library(ggplot2)
library(biomformat)
library(Biostrings)
library(ggpubr)
library(tictoc)

##### BasalMaxima analysis of primate data

B <- read_biom('data/primates/final.withtax.min10.biom')
Data <- biom_data(B)
taxonomy <- data.frame('sequence'=sapply(B$rows,getElement,'id'),
                       'taxonomy'=sapply(B$rows,FUN=function(x) paste(x$metadata$taxonomy,collapse='; ')))
taxonomy$taxonomy <- as.character(taxonomy$taxonomy)
X <- read.csv('data/primates/Reduced_mapping_4persp_nohyb_RC2036.txt',sep = '\t',stringsAsFactors = F)
names(X)[1] <- 'SampleID'
Data <- Data[,match(X$SampleID,colnames(Data))]
all.equal(colnames(Data),X$SampleID)



# trees -------------------------------------------------------------------
seqs <- as.list(rownames(Data))
paste(paste('>',seqs,sep=''),seqs,sep='\n') %>%
  write(file='data/primates/final_withtax_min10.fasta',sep='\n')
seqs <- Biostrings::DNAStringSetList(seqs)
write.dna(seqs,file='data/primates/final_withtax_min10.fasta')

col.tree <- read.tree('data/primates/furcated_host_tree.tre')
row.tree <- read.tree('data/primates/final_withtax_min10.tre')
all(row.tree$tip.label %in% rownames(Data))
all(col.tree$tip.label %in% colnames(Data))
Data <- Data[row.tree$tip.label,]
Data <- Data[,col.tree$tip.label]



# Is there signal? --------------------------------------------------------
W <- treeBasis(row.tree)
V <- treeBasis(col.tree)

U <- t(W) %*% Data %*% V

shuffleData <- function(Data,rows=TRUE,cols=TRUE){
  if (rows){
    ix.rows <- sample(nrow(Data))
  } else {
    ix.rows <- 1:nrow(Data)
  }
  if (cols){
    ix.cols <- sample(ncol(Data))
  } else {
    ix.cols <- 1:ncol(Data)
  }
  return(Data[ix.rows,ix.cols])
}

Unull <- t(W) %*% shuffleData(Data) %*% V


DF <- data.table('z'=c(log(c(as.matrix(U))^2),
                       log(c(as.matrix(Unull))^2)),
                 'Dataset'=rep(c('Primates','Null'),each=nrow(U)*ncol(U)))

ggplot(DF[z>-10],aes(z,by=Dataset,fill=Dataset))+
  geom_histogram(bins=100,alpha=0.4,position='identity')

z <- DF[Dataset=='Primates' & z>-10,z]
znull <- DF[Dataset=='Null' & z>-10,z]



# dendromap ---------------------------------------------------------------


tic()
dm <- dendromap(Data,row.tree,col.tree,Pval_threshold = 0.001)  
toc()

# par(mfrow=c(1,1))
# plot(ecdf(z),lwd=2)
# lines(ecdf(znull),lwd=2,
#       col=rgb(0.6,0.1,0.2,0.6))
# 
# z <- log(c(as.matrix(U))^2)
# z <- z[z>-10]
# x <- seq(min(z),max(z),length.out = 1e3)
# par(mfrow=c(1,3))
# hist(z,breaks=100,freq = F)
# lines(x,dnorm(x,mean=mean(z),sd=sd(z)),lwd=2,col='blue')
# 
# plot(ecdf(z))
# lines(ecdf(rnorm(1e4,mean=mean(z),sd=sd(z))),col='red',lwd=2)
# 
# hist(1-pnorm(z,mean=mean(z),sd=sd(z)),breaks=100)
# 
# 
# ### counts
# S <- basalMaxima(as.matrix(Data),row.tree,col.tree,ncores=7)
# 
# ### Presence-Absence
# Data[Data>0] <- 1
# V <- treeBasis(row.tree)
# W <- treeBasis(col.tree)
# U <- t(V) %*% Data %*% W
# 
# y <- c(as.matrix(U))
# z <- log(y^2)
# z <- z[z>-10]
# x <- seq(min(z),max(z),length.out = 1e3)
# par(mfrow=c(1,3))
# hist(z,breaks=100,freq = F)
# 
# S.PA <- basalMaxima(as.matrix(Data),row.tree,col.tree,threshold = 0,ncores=7)
# save(list=ls(),file='data/primates/basalMax_workspace')
# 
# # visualize ---------------------------------------------------------------
# 
# S[,orientation:=sign(value)]
#  
# ### code copied from pathViz
# colmap <- data.table('color'=viridis::viridis(length(unique(S$col.node))),
#                     'node'=unique(S$col.node))
# row.tree$edge.length[row.tree$edge.length>0.2]<-0.2
# rtr <- phytools::force.ultrametric(row.tree)
# par(mfrow=c(1,2))
# ape::plot.phylo(rtr,main='Microbial Tree, Count Analysis',
#                 show.tip.label = F,use.edge.length = T,edge.width = 0.01,edge.color = 'grey')
# nodelabels(text=rep(' ',nrow(S)),node=S$row.node,cex=0.5,
#           bg = colmap[match(S$col.node,node),color],frame = 'circle')
# 
# ape::plot.phylo(col.tree,main='Primate Tree',show.tip.label = F)
# nodelabels(text=rep(' ',nrow(colmap)),node=colmap$node,cex=1,bg=colmap$color,frame='circle')
# 
# row.depths <- node.depth.edgelength(row.tree)[S$row.node]
# col.tree$edge.length[is.na(col.tree$edge.length)] <- 0
# col.depths <- node.depth.edgelength(col.tree)[S$col.node]
# 
# plot(row.depths,col.depths)
# abline(coef(glm(col.depths~row.depths)))
# # Comparison to null simulations ------------------------------------------
# 
# rData <- function(m=1e3,n=10,mean=0,sd=1,mean_counts=1e4,presence.absence=F){
#   clr_inv <- function(x) exp(x)/sum(exp(x))
#   X <- matrix(rnorm(n=m*n,mean=mean,sd=sd),nrow=m) %>%
#     apply(2,clr_inv) %>% 
#     apply(2,FUN=function(p) rmultinom(1,rpois(1,lambda=mean_counts),p))
#   if (presence.absence){
#     X[X>0] <- 1
#   }
#   return(X)
# }
# Data <- biom_data(B)
# Data <- Data[row.tree$tip.label,col.tree$tip.label]
# PA.Data <- Data
# PA.Data[PA.Data>1] <- 1
# m <- nrow(Data)
# n <- ncol(Data)
# V <- treeBasis(row.tree)
# W <- treeBasis(col.tree)
# 
# mu=apply(Data,1,mean)
# mx <- max(mu)
# # X <- rData(m,n,mean=2*mu/mx,sd=3.5,presence.absence=F,mean_counts = mean(apply(Data,2,sum)))
# X <- rData(m,n,mean=log(mu),sd=1.5,presence.absence=F,mean_counts = mean(apply(Data,2,sum)))
# 
# 
# mns <- rowMeans(X)
# vs <- apply(X,1,var)
# 
# par(mfrow=c(1,1))
# plot(sort(mu/sum(mu),decreasing = T),xlab='Rank',ylab='Relative Abundance',log='y',type='l',lwd=2)
# lines(sort(mns/sum(mns),decreasing=T),lwd=2,col='blue')
# legend('topright',legend=c('Primate Data','Null Simulation'),lwd=c(2,2),col=c('black','blue'))
# 
# par(mfrow=c(1,2))
# plot(mu,apply(Data,1,var),xlab='mean',ylab='var',log='xy',main='Primate Data')
# abline(0,1)
# abline(1,2,col='blue',lwd=2)
# 
# plot(mns,vs,xlab='mean',ylab='var',log='xy',main='Null Data')
# abline(0,1)
# abline(1,2,col='blue',lwd=2)
# glm(log(vs)~log(mns)) ## slope >2
# 
# 
# plot(sort(rowSums(X),decreasing = T),log='y')
# 
# 
# png(filename = 'Figures/Primate_vs_Null_data_statistic_comparison_counts.png',height = 8,width=14,units='in',res = 200)
#   U <- t(V) %*% X %*% W
#   y <- c(U)
#   z <- log(y^2)[log(y^2)>-10]
#   x <- seq(min(z),max(z),length.out=1e3)
#   par(mfrow=c(2,2))
#   hist(z,breaks=100,freq = F,xlab='log(u^2)',main='Null Data',xlim=c(-10,15))
#   lines(x,dnorm(x,mean=mean(z),sd=sd(z)),col='blue')
#   qqplot(z,rnorm(1e3,mean=mean(z),sd=sd(z)),xlab='log(u^2)',ylab='Gaussian Quantiles')
#   abline(0,1)
#   
#   U <- t(V) %*% Data %*% W
#   y <- c(as.matrix(U))
#   z <- log(y^2)[log(y^2)>-10]
#   x <- seq(min(z),max(z),length.out=1e3)
#   
#   hist(z,breaks=100,freq = F,xlab='log(u^2)',main='Primate Data',xlim=c(-10,15))
#   lines(x,dnorm(x,mean=mean(z),sd=sd(z)),col='blue')
#   qqplot(z,rnorm(1e3,mean=mean(z),sd=sd(z)),xlab='log(u^2)',ylab='Gaussian Quantiles')
#   abline(0,1)
# dev.off()
# 
# 
# X.PA <- X
# X.PA[X.PA>0] <- 1
# png(filename = 'Figures/Primate_vs_Null_data_statistic_comparison_PA.png',height = 8,width=14,units='in',res = 200)
#   U <- t(V) %*% X.PA %*% W
#   y <- c(U)
#   z <- log(y^2)[log(y^2)>-10]
#   x <- seq(min(z),max(z),length.out=1e3)
#   par(mfrow=c(2,2))
#   hist(z,breaks=100,freq = F,xlab='log(u^2)',main='Null Data',xlim=c(-10,4))
#   lines(x,dnorm(x,mean=mean(z),sd=sd(z)),col='blue')
#   qqplot(z,rnorm(1e3,mean=mean(z),sd=sd(z)),xlab='log(u^2)',ylab='Gaussian Quantiles')
#   abline(0,1)
#   
#   U <- t(V) %*% PA.Data %*% W
#   y <- c(as.matrix(U))
#   z <- log(y^2)[log(y^2)>-10]
#   x <- seq(min(z),max(z),length.out=1e3)
#   
#   hist(z,breaks=100,freq = F,xlab='log(u^2)',main='Primate Data',xlim=c(-10,4))
#   lines(x,dnorm(x,mean=mean(z),sd=sd(z)),col='blue')
#   qqplot(z,rnorm(1e3,mean=mean(z),sd=sd(z)),xlab='log(u^2)',ylab='Gaussian Quantiles')
#   abline(0,1)
# dev.off()
# 
# 
# 
# U <- t(V) %*% X %*% W
# y <- c(U)
# z <- log(y^2)[log(y^2)>-10]
# DF <- data.frame("z"=z,"Dataset"='NULL')
# 
# U <- t(V) %*% Data %*% W
# y <- c(as.matrix(U))
# z <- log(y^2)[log(y^2)>-10]
# DF <- rbind(DF,data.frame("z"=z,"Dataset"='Primate'))
# 
# g1 <- ggplot(DF,aes(z,color=Dataset,fill=Dataset),alpha=0.5)+
#         geom_density(alpha=0.5)+
#         ggtitle('Statistic Distribution under Count Data')
# 
# U <- t(V) %*% X.PA %*% W
# y <- c(U)
# z <- log(y^2)[log(y^2)>-10]
# DF <- data.frame("z"=z,"Dataset"='NULL')
# 
# U <- t(V) %*% PA.Data %*% W
# y <- c(as.matrix(U))
# z <- log(y^2)[log(y^2)>-10]
# DF <- rbind(DF,data.frame("z"=z,"Dataset"='Primate'))
# 
# g2 <- ggplot(DF,aes(z,color=Dataset,fill=Dataset),alpha=0.5)+
#   geom_density(alpha=0.5)+
#   ggtitle('Statistic Distribution under Presence/Absence Data')
# 
# 
# ggarrange(g1,g2,ncol=2)
# ggsave('Figures/Primate_vs_Null_data_statistic_comparison.png')
# 
# 
# 
# # node depth distribution -------------------------------------------------
# nds <-Ntip(row.tree)+(1:Nnode(row.tree))
# DF <- data.frame('node'=c(rep('basalMax',times=nrow(S)),rep('tree',times=length(nds))),
#                  'depth'=c(node.depth.edgelength(row.tree)[S$row.node],
#                            node.depth.edgelength(row.tree)[nds]))
# 
# par(mfrow=c(1,1))
# plot(ecdf(DF$depth[DF$node=='tree']),xlab='Height',main='Node Height in tree vs. basalMax')
# lines(ecdf(DF$depth[DF$node=='basalMax']),col='blue',lwd=2)
# ks.test(DF$depth[DF$node=='tree'],DF$depth[DF$node!='tree'])
