### playing around with asr

rm(list=ls())
gc()
library(phylofactor)
library(dendromap)

# m=10
# n=5
# set.seed(1)
# row.tree <- rtree(m)
# col.tree <- rtree(n)
# 
# X <- rnorm(m*n) %>% matrix(nrow=m)
# rownames(X) <- row.tree$tip.label
# colnames(X) <- col.tree$tip.label
# 
# W <- treeBasis(row.tree)
# V <- treeBasis(col.tree)
# 
# Grps <- getPhyloGroups(row.tree)
# nds <- which(row.tree$edge[,2]>m)
# 
# edgeW <- getPhyloGroups(row.tree) %>% sapply(ilrvec,n=m)
# colnames(edgeW) <- row.tree$edge[,2]
# edgeV <- getPhyloGroups(col.tree) %>% sapply(ilrvec,n=n)
# 
# Crow <- vcv.phylo(row.tree)
# Ccol <- vcv.phylo(col.tree)
# 
# A <- apply(X,2,FUN=function(x,tr) phytools::fastAnc(tr,x),tr=row.tree)
# 
# a1=phytools::fastAnc(row.tree,X[,1])
# 
# ### tips 't10','t6','t9' , indexes 1:3, separate from the rest. This will be a contrast at node 11, contrasting nodes 12|14
# a1['12']
# a1['14']
# 
# 
# Aedge <- t(edgeW[,nds]) %*% X
# 
# plot(c(A[2:nrow(A),]),c(Aedge))
# 
# 
# 
# ### fastAnc uses ace(...,method='pic') with a tree rooted at the node
# par(mfrow=c(1,2))
# plot(row.tree)
# nodelabels()
# tt=multi2di(root(row.tree,node=18))
# plot(tt)
# nodelabels()
# 
# ace(X[,1],tt,method = 'pic')$ace[1]
# A['18',1]
# Aedge['18',1]
# 
# 
# a=ace(X[,1],tt,method = 'pic')
# a$ace[1]
# mean(a$CI95[1,])
# 
# ### we can use these to translate our ancestral character estimates to z-scores
# 
# ancz <- function(x,tree,node){
#   
#   a <- ape::ace(x,multi2di(root(tree,node=node)),method='pic')
#   (a$ace[1]-a$CI95[1,1])/qnorm(0.025)  ## under N(0,1), CI95 goes from
#   qnorm(0.025)
# }
# 



# bird biogeography -------------------------------------------------------

rm(list=ls())
gc()
load('data/birds/bird_dendromap_workspace')
m <- nrow(X)
n <- ncol(X)


W <- getPhyloGroups(row.tree) %>% sapply(ilrvec,n=nrow(X))
V <- getPhyloGroups(col.tree) %>% sapply(ilrvec,n=ncol(X))

# Cinv <- ape::vcv.phylo(row.tree) %>% chol %>% chol2inv
# Cinv <- chol2inv(chol(Crow))

# Xcor <- Cinv %*% X
# U <- t(W) %*% Xcor %*% V

U <- t(W) %*% X %*% V

## the biggest element here is 4083, 8; the next 4007, 8
ix=which(abs(U)>0) %>% arrayInd(.dim = c(nrow(U),ncol(U))) %>% as.data.table
names(ix) <- c('row.edge','col.edge')
ix$score <- U[abs(U)>0]
setkey(ix,row.edge,col.edge,score)

ix[,rank:=rank(-abs(score))]

setkey(ix,rank)

## comparison to null?
hist(ix$score,breaks=100)
Unull <- t(W) %*% X[sample(nrow(X)),sample(ncol(X))] %*% V
hist(c(Unull),breaks=100)

cdf <- ecdf(c(Unull))
Un2 <- t(W) %*% X[sample(nrow(X)),sample(ncol(X))] %*% V
hist(cdf(c(Un2)),breaks=10)
plot(ecdf(cdf(c(Un2))))

ix[,P:=1-cdf(score)]

setkey(ix,P)




##whittle down to one column assignment per edge - may even refine further based on score significance

gettips <- function(edge,tree)  phangorn::Descendants(tree,tree$edge[edge,2],'tips')[[1]] %>% tree$tip.label[.]

i=3
gettips(ix$row.edge[i],row.tree)
gettips(ix$col.edge[i],col.tree)


ix
gondwanas <- ix[col.edge==4 & score>0]
gettips(4083,row.tree)  #gondwn
gettips(4007,row.tree)

gettips(147,row.tree)
gettips(5,col.tree)

ix[row.edge==1]
##SA score 4.22
##GW score 2.016
#tinamous + rheas: row.edge 1(+) or 106 (-)
#While the S-American (col.edge 8) classification for tinamous + Rheas 
#dominates the Gondwanan classification (col.edge 4), both have P=0
# and one can hope the presence of an even stronger SAmerican descendant - row.edge 914
# could 

gettips(1,row.tree) #tinamous + rheas

row.tree$edge[1:3,]
gettips(2,row.tree) #tinamous - this is the S American one
ix[row.edge==2]  
#SA score 4.59 - one possibility is to use "basalMax" to consider a subset of edges
#GW score 1.81

### Uncorrected:
# row.edge col.edge      score  rank
# 1:        1        2 -1.9243051 32056
# 2:        1        1 -2.0163674 32664
# 3:        1        7 -1.4745173 32672
# 4:        1        3 -0.5452305 41368
# 5:        1        4  2.0163674 51056
# 6:        1        5 -0.2737973 51104
# 7:        1        6  2.2399220 51488
# 8:        1        8  4.2178503 51493

### Inverse-covariance corrected:
# row.edge col.edge         score  rank
# 1:        1        5 -0.0003622184 26801
# 2:        1        1 -0.0011636718 26802
# 3:        1        6  0.0014594218 26804
# 4:        1        7 -0.0004560421 26806
# 5:        1        3 -0.0006521043 26807
# 6:        1        4  0.0011636718 26810
# 7:        1        8  0.0022434615 26814
# 8:        1        2 -0.0007730968 26815

# filtering ix ------------------------------------------------------------





## only allow negative scores for basal edges - 1 and 4; also, just choose one. 
row.basal.edges <- which(row.tree$edge[,1]==(ape::Ntip(row.tree)+1))
col.basal.edges <- which(col.tree$edge[,1]==(ape::Ntip(col.tree)+1))

## allow pos/neg values of the first, and discard the second

ix[row.edge==row.basal.edges[1] | col.edge==col.basal.edges[2],score:=abs(score)]
ix <- ix[row.edge != row.basal.edges[2] & col.edge != col.basal.edges[2]]

ix <- ix[score>0 | col.edge==1]
ix <- ix[,list(col.edge=col.edge[which.max(abs(score))],
               score=score[which.max(abs(score))]),by=row.edge]
setkey(ix,score)

gondwanas <- ix[col.edge==1]$row.edge[1:10]  #negative scores for edge 1 are Gondwanan


## biggest gondwana clade
i=6
gettips(gondwanas[i])


# ix$score <- U[abs(U)>5]
# 
# [,1] [,2]
# [1,] 4007    8
# [2,] 4083    8
# [3,] 5171    8
# [4,] 6382    8
# [5,] 6388    8
# [6,] 6436    8

U[abs(U)>15]

# 22.82622  23.56856  17.63574 -15.39195 -15.36746 -15.17179

U[c(4083,4007),8]  #8 is South America - these are birds present in SA

## 4083 has the high score with col.edge 8, but 4007 is basal and includes the Pittas etc., which are also suboscines;
## 4083 and 4007 also have somewhat large scores with col.edge 6 - Af+SA.
## Which of these captures more variance in X?
edgeFstat <- function(X,lineage_table,W,V){
  w <- W[,lineage_table$row_edge,drop=F]
  v <- V[,lineage_table$col_edge,drop=F]
  Xhat <- w %*% (t(w) %*% X %*% v) %*% t(v)
  ix <- which(Xhat!=0) #indexes for blocks being predicted by wDv'
  n <- length(ix)
  rss <- sum((X-mean(X)-Xhat)[ix]^2)
  ess <- sum(c(Xhat[ix])^2)
  dfm <- ncol(w)+1
  df0 <- n-(ncol(w)+1)
  Fstat <- (ess/dfm)/(rss/df0)
  return(Fstat)
}
re=4083
ce=8
lineage_table=data.table('row_edge'=4007,
                         'col_edge'=8)
edgeFstat(X,lineage_table,W,V)

#4083,8   :   1988.012
#4083,6   :   529.2093
#4007,6   :   481.8902
# phangorn::Descendants(col.tree,col.tree$edge[gondwanas$col.edge[i],2],'tips')[[1]] %>% col.tree$tip.label[.]

s1=phangorn::Descendants(row.tree,row.tree$edge[4083,][2],'tips')[[1]] %>% row.tree$tip.label[.]
phangorn::Descendants(col.tree,col.tree$edge[8,][2],'tips')[[1]] %>% col.tree$tip.label[.]

s2=phangorn::Descendants(row.tree,row.tree$edge[4007,][2],'tips')[[1]] %>% row.tree$tip.label[.]

setdiff(s2,s1)


phangorn::Descendants(row.tree,row.tree$edge[5171,2],'tips')[[1]] %>% row.tree$tip.label[.]


#### functions needed:
#edgeTonode - convert edges in rc_table to nodes
#makeDescendantTable - w/o sign(stat)
#find_joinables - w/o sign(stat)



# Rheas -------------------------------------------------------------------

