####### seq search
library(igraph)
library(phylofactor)

set.seed(1)
m <- 200
n <- 10
row.tree <- rtree(m)
col.tree <- rtree(n)

W <- treeBasis(row.tree)
V <- treeBasis(col.tree)

### let's have a three node sequence contain an effect:
par(mfrow=c(1,2))
plot(row.tree)
nodelabels(c(340,341,344,353,350),c(340,341,344,353,350),bg='red')
plot(col.tree)
nodelabels(c(12,14,15,16,18),c(12,14,15,16,18),bg='red')

###
answer <- data.table('row.node'=c(340,341,344,353,350),
                     'col.node'=c(12,14,15,16,18),
                     'stat'=c(10,8,-9,9,11))

## orientations are: +,+,-,+,+
set.seed(2)
eta <- W[,answer$row.node-m] %*% diag(answer$stat) %*% t(V[,answer$col.node-n])
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


### make table for row-col node pairs
makeRCtable <- function(N,row.tree,col.tree,n_sim=NULL){
  W <- treeBasis(row.tree)
  V <- treeBasis(col.tree)
  U <- t(W) %*% N %*% V
  row.nodes <- sapply(rownames(U),strsplit,'_') %>% sapply(getElement,2) %>% as.numeric
  col.nodes <- sapply(colnames(U),strsplit,'_') %>% sapply(getElement,2) %>% as.numeric
  
  rc_table <- data.table(expand.grid('row.node'=row.nodes,
                                     'col.node'=col.nodes),
                         'stat'=c(U))
  rm('U')
  gc()
  if (is.null(n_sim)){
    n_sim <- ceiling(10000/(nrow(N)*ncol(N)))
  }
  rsim <- function(x,N,W,V) abs(c(t(W) %*% N[sample(nrow(N)),sample(ncol(N))] %*% V))
  null_cdf <- ecdf(sapply(1:n_sim,rsim,N,W,V))
  rc_table[,P:=1-null_cdf(abs(stat))]
  if (any(rc_table$P==0)){
    warning('Some P-values are 0. Will replace with 1/(n_sim*nrow(N)*ncol(N))')
    rc_table[P==0,P:=1/(2*n_sim*nrow(N)*ncol(N))]
  }
  rc_table[,rc_index:=1:.N]
  setkey(rc_table,row.node,col.node)
  return(rc_table)
}

rc_table <- makeRCtable(N,row.tree,col.tree)
# rc_table[order(P)][1:30]
# answer

### filter rc_table by P-val:
Pval_threshold=0.01
rc_table <- rc_table[P<=Pval_threshold]

## make tree of RC pair ancestor/descendant relations
## based on descendants of indexes

makeDescendantTable <- function(i=1,rc_table,
                                terminalRowNodes,
                                terminalColNodes,
                                Row_Descendants,
                                Col_Descendants){
  row.node <- rc_table$row.node[i]
  col.node <- rc_table$col.node[i]
  if (row.node %in% terminalRowNodes | 
      col.node %in% terminalColNodes){
    return(NULL)
  } else {
    sgn <- sign(rc_table$stat[i])
    if (sgn==1){ ## pos-pos and neg-neg
      rowdescPos <- Row_Descendants[[toString(row.node)]][['pos']]
      coldescPos <- Col_Descendants[[toString(col.node)]][['pos']]
      rowdescNeg <- Row_Descendants[[toString(row.node)]][['neg']]
      coldescNeg <- Col_Descendants[[toString(col.node)]][['neg']]
    } else {
      rowdescPos <- Row_Descendants[[toString(row.node)]][['pos']]
      coldescPos <- Col_Descendants[[toString(col.node)]][['neg']]
      rowdescNeg <- Row_Descendants[[toString(row.node)]][['neg']]
      coldescNeg <- Col_Descendants[[toString(col.node)]][['pos']]
    }
    ix <- rc_table[(row.node %in% rowdescPos & col.node %in% coldescPos) |
                     (row.node %in% rowdescNeg & col.node %in% coldescNeg)]$rc_index
    if (length(ix)>0){
      return(data.table('ancestor'=rc_table$rc_index[i],
                      'descendant'=ix))
    } else {
      return(NULL)
    }
  }
}

makeRCMap <- function(rc_table,row.nodemap,col.nodemap){
  row.nodes <- unique(rc_table$row.node)
  col.nodes <- unique(rc_table$col.node)
  Row_Descendants <- lapply(row.nodes,getIndexSets,row.nodemap) %>%
    lapply(FUN=function(x,a) lapply(x,intersect,a),a=row.nodes)
  Col_Descendants <- lapply(col.nodes,getIndexSets,col.nodemap) %>%
    lapply(FUN=function(x,a) lapply(x,intersect,a),a=col.nodes)
  
  names(Row_Descendants) <- row.nodes
  names(Col_Descendants) <- col.nodes
  
  terminalRowNodes <- row.nodes[sapply(Row_Descendants,FUN=function(x) length(unlist(x)))==0]
  terminalColNodes <- row.nodes[sapply(Row_Descendants,FUN=function(x) length(unlist(x)))==0]
  ## we don't have to look for descendants of these row/col nodes
  
  ix <- which((!rc_table$row.node %in% terminalRowNodes) &
              (!rc_table$col.node %in% terminalColNodes))
  maps <- lapply(ix,makeDescendantTable,
         rc_table,
         terminalRowNodes,
         terminalColNodes,
         Row_Descendants,
         Col_Descendants)
  RCmap <- rbindlist(maps)
  RCmap[,terminal:=!(descendant %in% ancestor)]
  setkey(RCmap,descendant,ancestor)
  return(RCmap)
}

RCmap <- makeRCMap(rc_table,row.nodemap,col.nodemap)


### The following function walks from leaves of the RCmap to 
### the most basal ancestors of those leaves within our RCmap.
rc_seqs <- function(i,RCmap){
  nd <- RCmap$descendant[i]
  
  ### note: a given node may have multiple paths towards the root
  ### e.g. it's possible to have nodes (1,2,3,4,nd)
  ### such that either (1,3,nd) or (2,4,nd) are sequences due to
  ### incompatible orientations of 1-2, 3-4, etc.
  ### Hence, for a given terminal node, we may need to make a list
  ### of possible sequences.
  ### junctures will be identifiable based on 4-nd and 3-nd relations
  ### but no 3-4 or 4-3 relation. At such a juncture, we'll replicate
  ### the rc.seq currently rooted at 'nd' and then walk up one juncture,
  ### and when we're done we'll return to walk up the other branch.
  ### We'll have to keep track of the index on node.seq and its un-walked branch on RCmap
  
    ancs <- sort(RCmap[descendant==nd,ancestor],decreasing = T)
    if (length(ancs)>1){
      ### Some of these ancs will be sequences
      ### others will be incompatible
      ### e.g. if we have 1-nd and 2-nd, 
      ### and 1-2 is in RCmap, then we know P(1,2,nd) dominates as a longer alignment
      ### thus we don't need to consider 1-nd branching.
      ### recall that lower indexes correspond to lower column nodes
      ### hence the highest index will not have descendants.
      jj=0
      seq <- NULL
      while(length(ancs)>0){
        searching=T
        jj=jj+1
        seq[[jj]] <- nd
        temp_ancs <- ancs
          while(searching==T){
            seq[jj] <- list(c(temp_ancs[1],seq[[jj]]))
            temp_ancs <- sort(RCmap[descendant==temp_ancs[1],ancestor],decreasing = T)
            if (length(temp_ancs)==0){
              searching=F
            }
          }
        ancs <- setdiff(ancs,seq[[jj]])
      }
    } else {
      seq <- list(c(ancs,nd))
    }
  return(seq)
}



check.joinable <- function(seq1,seq2,Seqs,rc_table,r.nodemap,c.nodemap){
  
  ### some of these lists are whole lineages
  ### some are alternative sequences of the same lineage
  
  ### e.g. 
  
  ## 1-2-5-6 and 1-2-3-4 are alternative seqs if they are not lineages
  ## they are lineages if the rc's of {3,5} don't share an 
  ## ancestor-descendant relation in the column/row tree.
  ## The best way to check this will be to find the row/col descendants of 3,5
  ## if their intersect is empty, they are compatible.
  
  ## Alternatively, if we have:
  ## 1-2-3-4 and 1-2-5-4, we know these are alternative seqs 
  ## since 3 and 5 are both on rootpath to 4.
  intrsct <- intersect(Seqs[[seq1]],Seqs[[seq2]])
  if (length(intrsct)==0){
    joinable <- FALSE
  } else {
    disagreements1 <- !Seqs[[seq1]]%in%intrsct
    disagreements2 <- !Seqs[[seq2]]%in%intrsct
    ### if we see FALSE TRUE FALSE - a true surrounded by falses, then these are alternatives
    if (all(c(1,-1) %in% diff(disagreements1)) | 
        all(c(1,-1) %in% diff(disagreements2))){ ### 
      joinable <- FALSE
    } else {
      ix1 <- Seqs[[seq1]][min(which(disagreements1))]
      ix2 <- Seqs[[seq2]][min(which(disagreements2))]
      
      
      nds1 <- unlist(rc_table[rc_index==ix1,c('row.node','col.node')])
      nds2 <- unlist(rc_table[rc_index==ix2,c('row.node','col.node')])
      if (nds1['row.node']==nds2['row.node']){
        joinable <- FALSE
      } else {
        RowDesc1 <- c(unlist(getIndexSets(nds1['row.node'],r.nodemap)),nds1['row.node'])
        RowDesc2 <- c(unlist(getIndexSets(nds2['row.node'],r.nodemap)),nds2['row.node'])
        Common_row_descendants <- intersect(RowDesc1,RowDesc2)
        ColDesc1 <- c(unlist(getIndexSets(nds1['col.node'],c.nodemap)),nds1['col.node'])
        ColDesc2 <- c(unlist(getIndexSets(nds2['col.node'],c.nodemap)),nds2['col.node'])
        Common_col_descendants <- intersect(ColDesc1,ColDesc2)
        n_common_row <- length(Common_row_descendants)
        n_common_col <- length(Common_col_descendants)
        if (n_common_row>0 & n_common_col>0){ ## descendants overlap ==> alternatives
          joinable <- FALSE
        } else { ## descendants don't overlap ==> compatible
          joinable <- TRUE
        }
      }
      
    }
  }
  return(joinable)
}


##### THE MAIN INTERNAL FUNCTION OF DENDROMAP
findLineages <- function(RCmap,rc_table,
                         row.nodemap.=row.nodemap,
                         col.nodemap.=col.noemap){
  ### will start with index as descendant and traverse up
  ix <- which(RCmap$terminal)
  Seqs <- lapply(ix,rc_seqs,RCmap) %>% unlist(recursive=FALSE) %>% unique
  n <- length(Seqs)
 
  
  tbl <- data.table('seq1'=rep(1:(n-1),times=(n-1):1),key='seq1')
  tbl[,seq2:=(seq1+1):n,by=seq1]
  
  joinability <- apply(t(tbl),2,FUN=function(x,s,rc,r,c) check.joinable(x[1],x[2],s,rc,r,c),
                         s=Seqs,rc=rc_table,r=row.nodemap,c=col.nodemap)
  tbl[,joinability:=joinability]
  ### we recursively group joinable sequences into unique index sets...
  ### is joinability transitive? No.
  ### However, it is true that the union of joinable sequences will dominate either sequence
  ### thus we can reduce the number of seqs we check grouping any mutually joinable sets.
  ### these sets will be complete subgraphs in the graph for joinability==TRUE
  ### We can find these with the function igraph::cliques
  joinables <- tbl[joinability==TRUE]
  joinable.edges <- split(joinables[,c(1,2)],seq(nrow(joinables))) %>% unlist
  G <- igraph::make_graph(joinable.edges,directed = F)
  joinable.seqs <- igraph::cliques(G) ### now we need to remove elements that are strict subsets of other sets
  joinable.seqs <- joinable.seqs[order(sapply(joinable.seqs,length),decreasing = F)]
  for (i in 1:(length(joinable.seqs)-1)){
    check.subset <- any(sapply(joinable.seqs[(i+1):length(joinable.seqs)],
                        FUN=function(a,b) all(b %in% a),b=joinable.seqs[[i]]))
    if (any(check.subset)){
      joinable.seqs[[i]] <- NA
    }
  }
  found.cliques <- !sapply(joinable.seqs,FUN=function(x) all(is.na(x)))
  joinable.seqs <- joinable.seqs[found.cliques]
  joinable.seqs <- c(setdiff(unique(unlist(tbl[,c(1,2)])),joinable.seqs),
                     joinable.seqs)
  get_rc_index <- function(x,Seqs.=Seqs) unique(unlist(Seqs[x]))
  Seqs <- lapply(joinable.seqs,get_rc_index)
  return(Seqs)
}

Lineages <- findLineages(RCmap,rc_table)

compute_score <- function(lineage,rc_table.=rc_table) rc_table[rc_index %in% lineage,-sum(log(P))]
scores <- sapply(Lineages,compute_score)

winner <- which.max(scores)

output_table <- rc_table[rc_index %in% Lineages[[winner]]]
answer
output_table


### now for each of these seqs, we get the rc indexes

