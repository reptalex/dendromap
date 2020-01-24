# dataset --------------------------------------------------------------

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
effect.size=5
noise_added=1
counts=1e4
compositional_link <- function(X,counts=1e4) apply(X,2,FUN=function(x,a) a*(exp(x)/sum(exp(x))),counts)
X <- (S$W %*% (effect.size*S$D) %*% t(S$V)) + rnorm(m*n,sd=noise_added)+baseline_abundance
N <- matrix(rpois(m*n,lambda=compositional_link(X,counts)),nrow=m,ncol=n)
# N <- matrix(rpois(m*n,lambda=exp(X)),nrow=m,ncol=n)
rownames(N) <- row.tree$tip.label
colnames(N) <- col.tree$tip.label

X <- N
# initialize --------------------------------------------------------------

# 
# X
# row.tree
# col.tree
ncores=NULL
# Pval_threshold=0.01        ## will tune min/max P-vals for F-stat search
maxP=0.01
nP=100
W=NULL
V=NULL
n_sim=NULL
nreps=20


# guts --------------------------------------------------------------------




base::cat(paste('Checking Data and tree compatibility'))
### Align dataset to trees
if (!all(rownames(X) %in% row.tree$tip.label)){
  stop('There are rownames(X) not in row.tree$tip.label')
} else {
  if (!all(row.tree$tip.label %in% rownames(X))){
    row.tree <- ape::drop.tip(row.tree,setdiff(rownames(X),row.tree$tip.label))
  }
  X <- X[row.tree$tip.label,]
}
if (!all(colnames(X) %in% col.tree$tip.label)){
  stop('There are colnames(X) not in col.tree$tip.label')
} else {
  if (!all(col.tree$tip.label %in% colnames(X))){
    col.tree <- ape::drop.tip(col.tree,setdiff(colnames(X),col.tree$tip.label))
  }
  X <- X[,col.tree$tip.label]
}
if (!is.null(ncores)){
  cl <- parallel::makeCluster(ncores)
  parallel::clusterEvalQ(cl,library(dendromap))
} else {
  cl <- NULL
}

base::cat(paste('\nMaking Nodemaps'))
row.nodemap <- dendromap:::makeNodeMap(row.tree)
col.nodemap <- dendromap:::makeNodeMap(col.tree)

base::cat(paste('\nMaking RC table with',n_sim,'null simulations'))
W <- treeBasis(row.tree)
V <- treeBasis(col.tree)
RC_table <- makeRCtable(X,row.tree,col.tree,W,V,n_sim)



# new functions -----------------------------------------------------------


# getBlocks <- function(w,v){
#   blocks <- vector('list',length=ncol(w))
#   for (i in 1:ncol(w)){
#     if (i==1){
#       
#     } else {
#       
#     }
#   } 
# }
rootpath <- function(node,tree){
  N=length(tree$tip.label)
  nds <- node
  i=1
  while(!(N+1) %in% nds){
    ix=tree$edge[,2]==nds[i]
    if (!any(ix)){
      stop(paste('Could not find find descendant node',nds[i],'in tree$edge'))
    } else {
      nds <- c(nds,tree$edge[ix,1])
      i=i+1
    }
  }
  return(nds)
}


findLineageTable <- function(Pval_threshold,rc_tbl,row.nodemap,col.nodemap,cl,nreps){
  if (nrow(rc_tbl)>1){
    row.nodes <- unique(rc_tbl$row.node)
    col.nodes <- unique(rc_tbl$col.node)
    Row_Descendants <- lapply(row.nodes,getIndexSets,row.nodemap) %>%
      lapply(FUN=function(x,a) lapply(x,intersect,a),a=row.nodes)
    Col_Descendants <- lapply(col.nodes,getIndexSets,col.nodemap) %>%
      lapply(FUN=function(x,a) lapply(x,intersect,a),a=col.nodes)
    
    names(Row_Descendants) <- row.nodes
    names(Col_Descendants) <- col.nodes
    RCmap <- tryCatch(makeRCMap(rc_tbl,Row_Descendants,Col_Descendants),error=function(e) NULL)
    if (is.null(RCmap) | all(RCmap$terminal)){
      rc_tbl[,Lineage:=1:.N]
      return(rc_tbl)
    } else {
      # base::cat(paste('\nRCmap has',nrow(RCmap),'rows'))
      Lineages <- find_lineages(RCmap,rc_tbl,Row_Descendants,Col_Descendants,cl,nreps)
      compute_score <- function(lineage,rc_tbl.=rc_tbl) rc_tbl[rc_index %in% lineage,-sum(log(P))]
      
      # base::cat(paste('\nRC-RC Tree traversal found',length(Lineages),'sequences of RCs. \n If this number is large, joining sequences by finding cliques will take a long time.'))
      # base::cat(paste('\nFiltering RC sequences into lineages'))
      
      i=0
      output <- NULL
      while (length(Lineages)>0){
        i=i+1
        
        scores <- sapply(Lineages,compute_score)
        winner <- which.max(scores)
        
        output_table <- rc_tbl[rc_index %in% Lineages[[winner]]]
        output_table[,Lineage:=i]
        output <- rbind(output,output_table)
        Lineages <- filter_winner(winner,Lineages,row.tree,rc_tbl,row.nodemap)
      }
      return(output)
    }
  } else {
    rc_tbl$Lineage <- 1
    return(rc_tbl)
  }
}

findOptimalP <- function(Pset,RC_table,row.nodemap,col.nodemap,cl,nreps){
  Fstats <- numeric(length(Pset))
  for (k in 1:length(Pset)){
    Pval_threshold <- Pset[k]
    rc_table <- RC_table[P<=Pval_threshold]
    lineage_table <- findLineageTable(Pval_threshold,rc_table,row.nodemap,col.nodemap,cl,nreps)
    Fstats[k] <- getFstat(X,lineage_table,W,V)
  }
  # plot(Pset,Fstats,type='l')
  k <- which.max(Fstats)
  return(Pset[k])
}

getFstat <- function(X,lineage_table,W,V){
  w <- W[,lineage_table$row.node-nrow(W),drop=F]
  v <- V[,lineage_table$col.node-nrow(V),drop=F]
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


# Fstat searching ---------------------------------------------------------


Pset <- unique(RC_table[P<=maxP,P]) %>% sort(decreasing=F)
if (length(Pset)>nP){
  Pset <- seq(min(Pset),maxP,length.out=nP)
}
Fstats <- numeric(length(Pset))
for (k in 1:length(Pset)){
  if (k==1){
    rc_table <- RC_table[P<=Pset[k]]
    lineage_table <- findLineageTable(Pset[k],rc_table,row.nodemap,col.nodemap,cl,nreps)
  } else {
    ## we only have to recompute lineages from rc_table if they are either descendants/ancestors/same-row.nodes
    ## of row.nodes in rc_table, OR if they are descendants of an ancestor of a new node (sister nodes).
    # We'll use rc_table of new row.nodes - rc_2 - to search for overlaps
    rc_2 <- RC_table[P>Pset[k-1] & P<=Pset[k]] 
    descendants <- lapply(rc_2$row.node,getIndexSets,row.nodemap) %>% unlist %>% unique
    ancestors <- sapply(rc_2$row.node,rootpath,row.tree) %>% unlist %>% c %>% unique
    ancestors_from_prev <- rc_table[row.node %in% ancestors,row.node]
    if (length(ancestors_from_prev)>0){
      cousins <- lapply(ancestors_from_prev,getIndexSets,row.nodemap) %>% unlist %>% unique
    } else {
      cousins <- NULL
    }
    affected_nds <- unique(c(ancestors,descendants,cousins))
    unaffected_lineages <- lineage_table[!row.node %in% affected_nds]
    rc_2 <- rbind(rc_2,rc_table[row.node %in% affected_nds])
    rc_2 <- rc_2[!duplicated(rc_2)]
    new_lineages <- findLineageTable(Pval_threshold,rc_2,row.nodemap,col.nodemap,cl,nreps)
    if (nrow(unaffected_lineages)>0){
      new_lineages$Lineage <- new_lineages$Lineage+max(unaffected_lineages$Lineage)
    }
    lineage_table <- rbind(unaffected_lineages,new_lineages)
    rc_table <- RC_table[P<=Pset[k]]
  }
  Fstats[k] <- getFstat(X,lineage_table,W,V)
}
plot(Pset,Fstats,type='l')
k <- which.max(Fstats)


OptimalP <- Pset[k]
# OptimalP <- findOptimalP(Pset,RC_table,row.nodemap,col.nodemap,cl,nreps)
output <- findLineageTable(OptimalP,RC_table[P<=OptimalP],row.nodemap,col.nodemap,cl,nreps)


# dk=3
# dum=findLineageTable(Pset[k+dk],RC_table[P<=Pset[k+dk]],row.nodemap,col.nodemap,cl,nreps)
# get winner --------------------------------------------------------------


# Pval_threshold <- Pset[k]
# rc_table <- RC_table[P<=Pval_threshold]
# if (nrow(rc_table)>1){
#   base::cat(paste('\n',nrow(rc_table),' RCs had P<=Pval_threshold at Pval_threshold=',Pval_threshold,sep=''))
#   
#   row.nodes <- unique(rc_table$row.node)
#   col.nodes <- unique(rc_table$col.node)
#   Row_Descendants <- lapply(row.nodes,getIndexSets,row.nodemap) %>%
#     lapply(FUN=function(x,a) lapply(x,intersect,a),a=row.nodes)
#   Col_Descendants <- lapply(col.nodes,getIndexSets,col.nodemap) %>%
#     lapply(FUN=function(x,a) lapply(x,intersect,a),a=col.nodes)
#   
#   names(Row_Descendants) <- row.nodes
#   names(Col_Descendants) <- col.nodes
#   RCmap <- tryCatch(makeRCMap(rc_table,Row_Descendants,Col_Descendants),error=function(e) NULL)
#   if (is.null(RCmap) | all(RCmap$terminal)){
#     output <- rc_table
#   } else {
#     base::cat(paste('\nRCmap has',nrow(RCmap),'rows'))
#     Lineages <- find_lineages(RCmap,rc_table,Row_Descendants,Col_Descendants,cl,nreps)
#     compute_score <- function(lineage,rc_table.=rc_table) rc_table[rc_index %in% lineage,-sum(log(P))]
#     
#     base::cat(paste('\nRC-RC Tree traversal found',length(Lineages),'sequences of RCs. \n If this number is large, joining sequences by finding cliques will take a long time.'))
#     base::cat(paste('\nFiltering RC sequences into lineages'))
#     
#     i=0
#     output <- NULL
#     while (length(Lineages)>0){
#       i=i+1
#       
#       scores <- sapply(Lineages,compute_score)
#       winner <- which.max(scores)
#       
#       output_table <- rc_table[rc_index %in% Lineages[[winner]]]
#       output_table[,Lineage:=i]
#       output <- rbind(output,output_table)
#       Lineages <- filter_winner(winner,Lineages,row.tree,rc_table,row.nodemap)
#     }
#   }
#   setkey(output,row.node,col.node)
# } else {
#   output <- rc_table
# }

if (!is.null(cl)){
  parallel::stopCluster(cl)
  rm('cl')
  gc()
}
output <- list('Lineages'=output,'Data'=X,'row.tree'=row.tree,'col.tree'=col.tree)
class(output) <- 'dendromap'