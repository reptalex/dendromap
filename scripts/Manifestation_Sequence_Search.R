### two-node searching
library(phylofactor)

set.seed(1)
m <- 30
n <- 10
row.tree <- rtree(m)
col.tree <- rtree(n)

W <- treeBasis(row.tree)
V <- treeBasis(col.tree)

### let's have a two node sequence contain an effect:
par(mfrow=c(1,2))
plot(row.tree)
nodelabels()
nodelabels(c(40,46),c(40,46),bg='red')
plot(col.tree)
nodelabels()
nodelabels(c(17,19),c(17,19),bg='red')

## column nodes 17, 19
## row nodes 40, 46
## orientation for (40,17): positive

set.seed(2)
eta <- W[,c(40,46)-m] %*% diag(c(10,-3)) %*% t(V[,c(17,19)-n])
par(mfrow=c(1,1))
image(t(eta))
X <- eta+matrix(rnorm(m*n),nrow=m)
rownames(X) <- row.tree$tip.label
colnames(X) <- col.tree$tip.label
clrinv <- function(x) exp(x)/sum(exp(x))
rmlt <- function(p,lambda=5e3) rmultinom(1,rpois(1,lambda),prob = p)
P <- apply(X,2,clrinv)
N <- apply(P,2,rmlt)
rownames(N) <- row.tree$tip.label
image(t(N))

# Node Maps ----------------------------------------------------------
#### make a node.map --> quick way to pull indexes for rows/columns of U corresponding to nodes' descendants

### compute the number of descendants on neg/pos side of each node
countDescendants <- function(node,tree,direction='pos'){
  children <- phangorn::Descendants(tree,node,type='children')
  D <- length(tree$tip.label)
  if (direction=='pos'){
    if (children[1]<=D){
      n <- 0
    } else {
      n <- length(setdiff(phangorn::Descendants(tree,children[1],type = 'all'),1:D))+1
    }
  } else {
    if (children[2] <=D){
      n <- 0
    } else {
      n <- length(setdiff(phangorn::Descendants(tree,children[2],type = 'all'),1:D))+1
    }
  }
  return(n)
}

makeNodeMap <- function(tree){
  D <- length(tree$tip.label)
  nodes <- (D+1):(D+ape::Nnode(tree))
  data.table('node'=nodes,
             'pos'=sapply(nodes,countDescendants,tree),
             'neg'=sapply(nodes,countDescendants,tree,direction='neg'),
             key='node') %>%
    return()
}
row.nodemap <- makeNodeMap(row.tree)
col.nodemap <- makeNodeMap(col.tree)



# Searching with Node Maps ------------------------------------------------

getIndexSets<- function(nd,nodemap){
  x <- vector(mode='list',length=2)
  x[[1]] <- nodemap[node==nd,(nd+1):(nd+pos)*(pos>0)]
  if (all(x[[1]]==0)){
    x[[1]] <- intersect(1,0)
  } 
  
  x[[2]] <- nodemap[node==nd,(nd+pos+1):(nd+pos+neg)*(neg>0)]
  if (all(x[[2]]==0)){
    x[[2]] <- intersect(1,0)
  }
  names(x) <- c('pos','neg')
  return(x)
}

U <- t(W) %*% N %*% V

## since we're looking for node sequences, we can focus on 
## 1) nodes that have descendants, i.e. nodemap[pos+neg>0]
## 2) for node pair with +orientation, must have 
#      row.nodemap[node==nd,c('pos','neg')]*col.nodemap[node==col.nd,c('pos','neg')]>0
#      row.nodemap[node==nd,c('pos','neg')]*col.nodemap[node==col.nd,c('neg','pos')]>0 for -orientation
#     This criterion ensures we're not searching for pos-descendants on row.tree if there are no pos.descendants on col.tree.

findBestSeq <- function(RowCol,U,row.tree,col.tree,row.nodemap,col.nodemap){
  row.node <- as.numeric(strsplit(rownames(U)[RowCol[1]],'_')[[1]][2])
  col.node <- as.numeric(strsplit(colnames(U)[RowCol[2]],'_')[[1]][2])
  m <- length(row.tree$tip.label)
  n <- length(col.tree$tip.label)
  z <- U[row.node-m,col.node-n]
  orientation <- c('neg','pos')[as.numeric(z>0)+1]
  
  row.nds <- getIndexSets(row.node,row.nodemap)
  col.nds <- getIndexSets(col.node,col.nodemap)
  if (orientation=='pos'){
    chk <- row.nodemap[node==row.node,c('pos','neg')]*col.nodemap[node==col.node,c('pos','neg')]
    if (all(chk==0)){
      M <- NULL
    } else if (all(chk>0)){
      M1 <- U[row.nds[['pos']]-m,col.nds[['pos']]-n,drop=F]
      M2 <- U[row.nds[['neg']]-m,col.nds[['neg']]-n,drop=F]
      ## change the following if statement & omega calculation if inputting P-value matrix
      ## lines 118-142
      m1 <- max(abs(M1))
      m2 <- max(abs(M2))
      if (m1>m2){
        M <- M1
      } else {
        M <- M2
      }
      rm(list=c('M1','M2','m1','m2'))
    } else if (chk$pos>0){
      M <- U[row.nds[['pos']]-m,col.nds[['pos']]-n,drop=F]
    } else {
      M <- U[row.nds[['neg']]-m,col.nds[['neg']]-n,drop=F]
    }
  } else {
    chk <- row.nodemap[node==row.node,c('pos','neg')]*col.nodemap[node==col.node,c('neg','pos')]
    if (all(chk==0)){
      M <- NULL
    } else if (all(chk>0)){
      M1 <- U[row.nds[['pos']]-m,col.nds[['neg']]-n,drop=F]
      M2 <- U[row.nds[['neg']]-m,col.nds[['pos']]-n,drop=F]
      ## change the following if statement & omega calculation if inputting P-value matrix
      ## lines 118-142
      m1 <- max(abs(M1))
      m2 <- max(abs(M2))
      if (m1>m2){
        M <- M1
      } else {
        M <- M2
      }
      rm(list=c('M1','M2','m1','m2'))
    } else if (chk$pos>0){
      M <- U[row.nds[['pos']]-m,col.nds[['neg']]-n,drop=F]
    } else {
      M <- U[row.nds[['neg']]-m,col.nds[['pos']]-n,drop=F]
    }
  }
  if (is.null(M) | any(dim(M)==0)){
    Seq <- NULL
  } else {
    ix <- which.max(abs(M))
    i <- ix %% nrow(M)
    if (i==0){
      i <- nrow(M)
    }
    row.node2 <- as.numeric(strsplit(rownames(M)[i],'_')[[1]][2])
    col.node2 <- as.numeric(strsplit(colnames(M)[ceiling(ix/nrow(M))],'_')[[1]][2])
    z2 <- M[i,ceiling(ix/nrow(M))]
    Seq <- data.table('row.node'=c(row.node,row.node2),'col.node'=c(col.node,col.node2),
                      'statistic'=c(z,z2))
  }
  return(Seq)
}

## e.g. the nodes below will be removed if we find a 40-46 sequence in row tree
forbidden.nodes <- nodepath(row.tree,40,m+1)
makeRowColSet <- function(U,forbidden.nodes=NULL){
  nds <- as.numeric(sapply(rownames(U),FUN=function(x) strsplit(x,'_')[[1]][2]))
  rows <- setdiff(1:nrow(U),which(nds %in% forbidden.nodes))
  RCset <- expand.grid(rows,1:ncol(U))
  RCset <- split(x = RCset,seq(nrow(RCset))) %>% lapply(unlist) %>%
    return()
}
## once we find a lineage, we'll disallow its ancestors on the microbial tree - this will change RowColSet to setdiff(1:nrow(U),forbidden.rows)

forbidden.nodes <- row.nodemap[pos+neg==0]$node
RowColSet <- makeRowColSet(U,forbidden.nodes)

seqs <- lapply(RowColSet,findBestSeq,U=U,
                row.tree=row.tree,col.tree=col.tree,
                row.nodemap=row.nodemap,col.nodemap=col.nodemap)

seqstat <- function(Seq){
  if (is.null(Seq)){
    return(-Inf)
  } else {
    return(sum(abs(Seq$statistic)))
  }
}

which.max(sapply(seqs,seqstat))

seqs[[122]]





Sim <- basalMaxima(N,row.tree,col.tree,W,V,2)

### Now we have two nodes, an ancestor & its descendant. We need to fill out
### the ancestors & all other descendants
max.nodeseq <- which.max(sapply(seqs,seqstat))

####### find ancestors
row.node <- seqs[[max.nodeseq]]$row.node[1]
col.node <- seqs[[max.nodeseq]]$col.node[1]

ancestors <- matchAncestors(U,row.node,col.node,row.tree,col.tree)


