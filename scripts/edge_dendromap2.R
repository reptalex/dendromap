rm(list=ls())
gc()


library(phylofactor)
# library(dendromap)
library(parallel)
load('data/birds/bird_dendromap_workspace')
source('R/edge_dendromap_fcns.R')
m <- nrow(X)
n <- ncol(X)



# edge rc_table -----------------------------------------------------------
# max_Pval=0.0001
# col_shuffle_only=TRUE
# W <- getPhyloGroups(row.tree) %>% sapply(ilrvec,n=nrow(X))
# V <- getPhyloGroups(col.tree) %>% sapply(ilrvec,n=ncol(X))
# 
# # Cinv <- ape::vcv.phylo(row.tree) %>% chol %>% chol2inv
# # Cinv <- chol2inv(chol(Crow))
# # Xcor <- Cinv %*% X
# # U <- t(W) %*% Xcor %*% V
# 
# U <- t(W) %*% X %*% V
# if (col_shuffle_only){
#   Unull <- t(W) %*% X[,sample(n)] %*% V
# } else {
#   Unull <- t(W) %*% X[sample(m),sample(n)] %*% V
# }
# 
# cdf <- ecdf(c(Unull))
# 
# 
# ix=which(abs(U)>0) %>% arrayInd(.dim = c(nrow(U),ncol(U))) %>% as.data.table
# names(ix) <- c('row.edge','col.edge')
# ix$stat <- U[abs(U)>0]
# setkey(ix,row.edge,col.edge,stat)
# 
# ix[,rank:=rank(-abs(stat))]
# ix[,P:=1-cdf(stat)]
# 
# ### need to modify basal edges
# # row.basal.edges <- which(row.tree$edge[,1]==(ape::Ntip(row.tree)+1))
# # col.basal.edges <- which(col.tree$edge[,1]==(ape::Ntip(col.tree)+1))
# 
# ## allow pos/neg values of the first, and discard the second
# 
# # ix[row.edge==row.basal.edges[1] | col.edge==col.basal.edges[1],stat:=abs(stat)]
# # ix <- ix[row.edge != row.basal.edges[2] & col.edge != col.basal.edges[2]]
# 
# ix <- ix[stat>0]
# 
# 
# setkey(ix,stat)
# 
# 
# rc_table <- ix
# 
# y <- log(c(Unull[Unull>0]))
# if (any(rc_table$P==0)){
#   nn=sum(rc_table$P==0)
#   base::cat(paste('\n',nn,' P-values were 0. Will estimate tail probabilities assuming log(stat^2)~rnorm(mu,sd)',sep=''))
#   y <- y[y>-20]
#   mu <- mean(y)
#   sig <- sd(y)
#   min.P <- rc_table[P>0,min(P)]
#   stat.min <- max(log(rc_table[P>0][P==min(P),stat]))
#   estimate_tail <- function(y1,mu,sig,ymin,pmin) (1-pnorm(y1,mu,sig))/(1-pnorm(ymin,mu,sig))*pmin
#   rc_table[P==0,P:=estimate_tail(log(stat),mu,sig,stat.min,min.P)]
# }
# 
# rc_table <- rc_table[P<=max_Pval]
# rc_table[,rc_index:=1:.N]


rc_table <- edge_rc_table(X,row.tree,col.tree,maxPval=max_Pval)

# edgeMap for quick descendant calculation -------------------------------------------------------------------

### note: we're now including tips of our tree (yay!)
### we'll need to consider this when using row.nodemap.
# makeNodeMap <- function(tree){
#   D <- length(tree$tip.label)
#   nodes <- (D+1):(D+ape::Nnode(tree))
#   data.table('node'=nodes,
#              'pos'=sapply(nodes,countDescendants,tree),
#              'neg'=sapply(nodes,countDescendants,tree,direction='neg'),
#              key='node') %>%
#     return()
# }

set.seed(1)
tree <- rtree(20)
plot(tree)
edgelabels()
makeEdgeMap <- function(tree){
  ## while nodemap had pos/neg, this will have edge,n = number descendants
  ntip <- ape::Ntip(tree)
  dt <- data.table('edge'=1:nrow(tree$edge),'n'=0)
  setkey(dt,edge)
  
  root_edges <- which(tree$edge[,1]==(ape::Ntip(tree)+1))
  
  tips <- data.table('edge'=which(tree$edge[,2]<=length(tree$tip.label)))
  tips[,anc:=match(tree$edge[edge,1],tree$edge[,2])]
  ## we'll loop through tips, joining them, defining their ancestor as tip.
  while(nrow(tips)>0){
    common_ancs <- as.numeric(names(which(table(tips$anc)==2)))
    jt <- tips[anc %in% common_ancs,list(edge=unique(anc),
                                         d1=edge[1],
                                         d2=edge[2]),by=anc]
    
    dt[match(jt$edge,edge),n:=n+dt[match(jt$d1,edge),n]+dt[match(jt$d2,edge),n]+2]
    tips <- tips[! anc %in% common_ancs]
    new_tips <- data.table('edge'=setdiff(common_ancs,root_edges))
    new_tips[,anc:=match(tree$edge[edge,1],tree$edge[,2])]
    tips <- rbind(tips,new_tips)
  }
  setkey(dt,edge)
  return(dt)
}
  
rowEdgeMap <- makeEdgeMap(row.tree)
colEdgeMap <- makeEdgeMap(col.tree)

# RC_table composition ----------------------------------------------------

row_edges = table(rc_table$row.edge) %>% sort(decreasing=T)
length(row_edges) ## 2660 edges in total
sum(row_edges>1) ## 732 edges appear more than once
sum(row_edges[row_edges>1]-1) ## 9945 of the 2660 edges arise from 2-4x counting of the same row edge
### Should we use a given row edge only once?
### What algorithm can whittle down a most-likely-valuable alignment?
get_genus <- function(spp,uniques=T){
  genera=strsplit(spp,'_') %>% sapply(getElement,1)
  if (uniques){
    return(unique(genera))
  } else{
    return(genera)
  }
}  

genera_to_edge <- function(gns,row.tree.=row.tree){
  spp=row.tree$tip.label
  genera <- get_genus(spp,uniques=F)
  tips <- match(gns,genera) %>% unique
  nd_mrca <- ape::getMRCA(row.tree,row.tree$tip.label[tips])
  edge <- which(row.tree$edge[,2]==nd_mrca)
  return(edge)
}


############ CASE STUDY: RHEAS ###############
## Rheas: edge 1 (Rhea) and 106 (non-Rhea) Need: 10
rc_table[row.edge==1]
# row.edge col.edge     stat   rank            P rc_index
# 1:        1        4 2.016367 3409.5 5.911509e-06     1486 #Gondwana
# 2:        1        6 2.239922 2681.5 4.650639e-06     1935 #SA/Africa
# 3:        1        8 4.217850  615.5 8.830377e-07     3319 #SA

### SA/Africa is more significant because only Emus/Cassowaries live on Australia. 
### We'd pick Gondwana if and only if we ID BOTH an Africa AND 

rheas <- row.tree$tip.label[phangorn::Descendants(row.tree,8106,'tips')[[1]]]
rhea.tree <- drop.tip(row.tree,setdiff(row.tree$tip.label,rheas))

par(mfrow=c(1,2))
plot(rhea.tree)
edgelabels(1+1:(ape::Nedge(rhea.tree)))
plot(col.tree)
edgelabels()

rc_table[row.edge %in% rowEdgeMap[edge==1,seq(edge,edge+n)]]
rc_table[row.edge==1]

rc_table[row.edge==1+95] ## cassowaries & Emus - Australia
rc_table[row.edge==1+1]  ## tinamous - south america
rc_table[row.edge==1+94]  ## tinamous - south america

focal_edge=1
descendants <- rowEdgeMap[edge==focal_edge,setdiff(seq(edge,edge+n),edge)]
rhea_table <- rc_table[row.edge %in% c(descendants,focal_edge)]


# RC_map ------------------------------------------------------------------


#' make rc-rc mapping for ancestor-descendant relationships & rc tree traversal
#' @export
#' @param rc_table made from \code{\link{makeRCtable}}
#' @param Row_Descendants named \code{getIndexSets} of all row.nodes
#' @param Col_Descendants named \code{getIndexSets} of all col.nodes
# makeRCMap <- function(rc_table,Row_Descendants,Col_Descendants){
#   row.nodes <- as.numeric(names(Row_Descendants))
#   col.nodes <- as.numeric(names(Col_Descendants))
#   
#   terminalRowNodes <- row.nodes[sapply(Row_Descendants,FUN=function(x) length(unlist(x)))==0]
#   terminalColNodes <- row.nodes[sapply(Row_Descendants,FUN=function(x) length(unlist(x)))==0]
#   ## we don't have to look for descendants of these row/col nodes
#   
#   ix <- which((!rc_table$row.node %in% terminalRowNodes) &
#                 (!rc_table$col.node %in% terminalColNodes))
#   maps <- lapply(ix,makeDescendantTable,
#                  rc_table,
#                  terminalRowNodes,
#                  terminalColNodes,
#                  Row_Descendants,
#                  Col_Descendants)
#   RCmap <- rbindlist(maps)
#   RCmap[,terminal:=!(descendant %in% ancestor)]
#   setkey(RCmap,descendant,ancestor)
#   return(RCmap)
# }
row.edges <- 1:ncol(W)
col.edges <- 1:ncol(V)
Row_Descendants <- rowEdgeMap[,list(list(setdiff(seq(edge,edge+n),edge))),by=edge]$V1
Col_Descendants <- colEdgeMap[,list(list(setdiff(seq(edge,edge+n),edge))),by=edge]$V1
makeDescendantTable <- function(i=1,rc_table,
                                terminalRowNodes,
                                terminalColNodes,
                                Row_Descendants,
                                Col_Descendants,
                                method='node'){
  if (method=='node'){
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
      ixs <- rc_table[(row.node %in% rowdescPos & col.node %in% coldescPos) |
                        (row.node %in% rowdescNeg & col.node %in% coldescNeg)]$rc_index
      if (length(ix)>0){
        return(data.table('ancestor'=rc_table$rc_index[i],
                          'descendant'=ixs))
      } else {
        return(NULL)
      }
    }
  } else {
    row.edge <- rc_table$row.edge[i]
    col.edge <- rc_table$col.edge[i]
    if (row.edge %in% terminalRowNodes|
        col.edge %in% terminalColNodes){
      return(NULL)
    } else {
      
      rowDesc <- Row_Descendants[[row.edge]]
      colDesc <- Col_Descendants[[col.edge]]
      
      ixs <- rc_table[(row.edge %in% rowDesc) & (col.edge %in% colDesc)]$rc_index
      if (length(ixs)>0){
        return(data.table('ancestor'=rc_table[i]$rc_index,
                          'descendant'=ixs))
      } else {
        return(NULL)
      }
    }
  }
}
edgeRCMap <- function(rc_table,Row_Descendants,Col_Descendants){

  row.terminals <- which(sapply(Row_Descendants,FUN=function(x) length(x)==0))
  col.terminals <- which(sapply(Col_Descendants,FUN=function(x) length(x)==0))

  ix <- which((!rc_table$row.edge %in% row.terminals) &
                (!rc_table$col.edge %in% col.terminals))
  maps <- lapply(ix,makeDescendantTable,
                 rc_table,
                 row.terminals,
                 col.terminals,
                 Row_Descendants,
                 Col_Descendants,
                 method='edge')
  RCmap <- rbindlist(maps)
  RCmap[,terminal:=!(descendant %in% ancestor)]
  setkey(RCmap,descendant,ancestor)
  return(RCmap)
}

## for testing makeDescendantTable
# for (i in 1:length(ix)){
#   makeDescendantTable(ix[i],rc_table,row.terminals,col.terminals,Row_Descendants,Col_Descendants,method='edge')
# }

RCmap <- edgeRCMap(rc_table,Row_Descendants,Col_Descendants)

ancs <- unique(RCmap$ancestor)
desc <- unique(RCmap$descendant)
basal_edges <- setdiff(ancs,desc)

# Rhea analysis -----------------------------------------------------------


# ancs <- unique(rhea_map$ancestor)
# desc <- unique(rhea_map$descendant)
# basal_edges <- setdiff(ancs,desc)

# rhea_map <- edgeRCMap(rhea_table,Row_Descendants,Col_Descendants)
# rhea_nds <- unique(rhea_map[terminal==TRUE,descendant])
# rhea_sqs <- lapply(rhea_nds,rc_seqs,rhea_map) %>% unlist(recursive=FALSE) %>% unique
# 
# 
# rhea_joinables <- find_joinables(rhea_sqs,rhea_table,Row_Descendants,Col_Descendants)

edge_ancestors <- function(edge,tree){
  rt=ape::Ntip(tree)+1
  anc=tree$edge[edge,1]
  if (anc==rt){
    return(NULL)
  } else {
    anc_edges <- NULL
    while(!anc==rt){
      new_anc_edge <- which(tree$edge[,2]==anc)
      anc_edges <- c(anc_edges,new_anc_edge)
      anc=tree$edge[new_anc_edge,1]
    }
    return(anc_edges)
  }
}

######## option 1: greedyColTree
incompatible_descendants <- function(rc_ix,rc_tbl,Row_Descendants,Col_Descendants,row.tree.=row.tree,col.tree.=col.tree){
  row.edg <- rc_tbl[rc_index==rc_ix,row.edge]
  col.edg <- rc_tbl[rc_index==rc_ix,col.edge]
  ### incompatible lineages will have either
  ## the same (or ancestral) col.edges in descendant row.edges
  ## or ancestral row.edges
  
  ## in other words, row_ancs are excluded and only descendant col.edges are allowed in descednant row edges
  
  row_descs <- Row_Descendants[[row.edg]]
  col_descs <- Col_Descendants[[col.edg]]
  row_ancs <- c(edge_ancestors(row.edg,row.tree),row.edg)
  incompatibles <- rc_tbl[row.edge %in% row_ancs,rc_index]
  if (is.null(col_descs)){ ## No possible descendants
    incompatibles <- c(incompatibles,rc_tbl[row.edge %in% row_descs,rc_index])
  } else {
    incompatibles <- c(incompatibles,rc_tbl[row.edge %in% row_descs & !col.edge %in% col_descs,rc_index])
  }
  return(incompatibles)
}



clean_sisters <- function(lineage,row.tree,rc_table){
  ### any sister row.edges with same col.edge will reclassify their mother as their col.edge
  repeated_col_edges <- as.numeric(names(which(table(lineage$col.edge)>1)))
  if (length(repeated_col_edges)==0){
    return(lineage)
  } else {
    row_edges <- lineage[col.edge %in% repeated_col_edges,row.edge]
    ancs <- row.tree$edge[row_edges,1]
    while(any(table(ancs)>1)){
      common_ancs <- as.numeric(names(which(table(ancs)==2)))
      for (anc in common_ancs){
        anc_edge <- which(row.tree$edge[,2]==anc)
        desc_edges <- row_edges[ancs==anc]
        col_edge <- lineage[row.edge==desc_edges[1],col.edge]
        lineage <- lineage[! row.edge %in% c(anc_edge,desc_edges)]
        if (nrow(rc_table[row.edge==anc_edge & col.edge==col_edge])==1){
          lineage <- rbind(lineage,rc_table[row.edge==anc_edge & col.edge==col_edge])
        } else {
          lineage <- rbind(lineage,
                           data.table('row.edge'=anc_edge,'col.edge'=col_edge,
                                      'stat'=NA,'rank'=NA,'P'=NA,'rc_index'=NA))
        }
      }
      row_edges <- lineage[col.edge %in% repeated_col_edges,row.edge]
      ancs <- row.tree$edge[row_edges,1]
    }
    return(lineage)
  }
}

filter_rc_tbl <- function(rc_tbl,Row_descendants.=RowDescendants){
  ### Filtering ### 
  ## We'd like to do this as few times as possible.
  ## One option is to start with basal indexes of rc_table for the min P-value cutoff
  ## and then use that filtered rc_table's basal-index shards for greedyColTree input
  ## When adding more rc-indexes with higher P-values, 
  ## we can simply check whether/not they, too, need to be excluded.
  
  ## RULES OF FILTERING:
  ## remove col.edges with more significant values among descendants (implies signal lost with inclusion of sister taxa)
  ### Implications:
  ## addition of P-values with higher values will not upset a lower P-value's filtering. 
  setkey(rc_tbl,row.edge,P)
  row.edges <- unique(rc_tbl$row.edge)

  for (ee in row.edges){
    dum <- rc_table[row.edge==ee,c('col.edge','P','rc_index')]
    if (nrow(dum)>0 & length(Row_Descendants[[ee]])>0){
      descs_tbl <- rc_table[row.edge %in% Row_Descendants[[ee]] & col.edge %in% dum$col.edge]
      descs <- descs_tbl[,list(Pmin=min(P)),by=col.edge]
      setkey(dum,col.edge)
      setkey(descs,col.edge) ## this may not contain all col.edges in dum
      dd <- descs[dum]
      removers <- dd[P>Pmin,rc_index]
      rc_tbl <- rc_tbl[!rc_index %in% removers]
    }
  }
  return(rc_tbl)
}



greedyColTree <- function(rc_tbl,RCmap.=RCmap,rc_table.=rc_table,row.tree.=row.tree,Row_Descendants.=Row_Descendants,filter=TRUE){
  ### RULES:
  ## 1) remove col.edges with more significant values among descendants (implies signal lost with inclusion of sister taxa)
  ##     see function filter_rc_tbl
  ## 2) pick most basal row.edge, choosing col.edge with lowest P-value  -- from filtering, this is guaranteed to be greater than that of any descendant for same col.edge
  ## 3) remove row.edge and all descendant rc_indexes with same col.edge from rc_tbl
  ## 4) repeat 2-3 until rc_tbl is empty
  
  if (filter){
    rc_tbl <- filter_rc_tbl(rc_tbl)
  }
  
 
  lineage=NULL
  done=F
  n=0
  setkey(rc_tbl,row.edge,P)
  while(!done){
    n=n+1
    lineage <- rbind(lineage,rc_tbl[1,]) ## next most basal row edge, and its lowest P-value col.edge
    ix <- rc_tbl$rc_index[1]
    incompatibles <- c(ix,incompatible_descendants(ix,rc_tbl,Row_Descendants,Col_Descendants))
    rc_tbl <- rc_tbl[!rc_index %in% incompatibles]
    if (nrow(rc_tbl)==0){
      done=TRUE
    }
  }
  lineage <- clean_sisters(lineage,row.tree,rc_table)
  lineage[,lineage_id:=rc_index[1]]
  return(lineage)
}

lineage <- greedyColTree(rhea_table,RCmap,rc_table,row.tree)



getLineage <- function(basal_ix,RCmap.=RCmap,rc_table.=rc_table,row.tree.=row.tree,
                       Row_Descendants.=Row_Descendants,Col_Descendants.=Col_Descendants){
  basal_row_edge <- rc_table[rc_index==basal_ix,row.edge]
  basal_col_edge <- rc_table[rc_index==basal_ix,col.edge]
  row_descs <- Row_Descendants[[basal_row_edge]]
  col_descs <- Col_Descendants[[basal_col_edge]]
  rc_tbl <- rc_table[row.edge %in% c(basal_row_edge,row_descs) & 
                       col.edge %in% c(basal_col_edge,col_descs)]
  lineage <- greedyColTree(rc_tbl,RCmap,rc_table,row.tree)
  lineage[,lineage_id:=basal_ix]
  setkey(lineage,col.edge,P)
  return(lineage)
}

### NOTE: I can filter rc_table first - for every basal_ix, run filter_rc_table then rbindlist
Simplify_rc_table <- function(rc_table,RCmap.=RCmap,
                             Row_Descendants.=Row_Descendants,
                             Col_Descendants.=Col_Descendants,
                             cl=NULL){
  ancs <- unique(RCmap$ancestor)
  descs <- unique(RCmap$descendant)
  basal_indexes <- setdiff(ancs,descs) ## ancestors w/o descendant
  
  ### We'll filter rc_table to: (1) only edges with ancestor/descendant relation
  ### i.e. must be in RCmap
  rc_table <- rc_table[rc_index %in% c(ancs,descs)]
  setkey(rc_table,row.edge,P)
  rc_table <- filter_rc_tbl(rc_table)
  return(rc_table)
}

rc_table <- Simplify_rc_table(rc_table)
RCmap <- edgeRCMap(rc_table,Row_Descendants,Col_Descendants)
# choosing edges to input to greedyColTree -----------------------------------------

# getBasalEdges <- function(rc_tbl,Row_Descendants){
#   setkey(rc_table,row.edge)
#   basal_edges <- rc_tbl[1,row.edge]
#   rc_tbl <- rc_tbl[!row.edge %in% c(basal_edges[1],Row_Descendants[[basal_edges[1]]])]
#   n=1
#   while(nrow(rc_tbl)>0){
#     n=n+1
#     basal_edges <- c(basal_edges,rc_tbl[1,row.edge])
#     rc_tbl <- rc_tbl[!row.edge %in% c(basal_edges[n],Row_Descendants[[basal_edges[n]]])]
#   }
#   return(basal_edges)
# }
# 
# basal_edges <- getBasalEdges(rc_table,Row_Descendants)
### problem: grabs two root edges - rheas and the other sister clade at root of bird tree...

ancs <- unique(RCmap$ancestor)
descs <- unique(RCmap$descendant)

basal_indexes <- setdiff(ancs,descs)

# lineages <- lapply(basal_indexes,getLineage,
#                         RCmap,rc_table,row.tree,Row_Descendants)

cl <- parallel::makeCluster(7)

parallel::clusterExport(cl,varlist=c('incompatible_descendants','clean_sisters','greedyColTree','getLineage','edge_ancestors',
                                     'RCmap','rc_table','Row_Descendants','Col_Descendants','row.tree','col.tree'))
parallel::clusterEvalQ(cl,library(data.table))
parallel::clusterEvalQ(cl,library(magrittr))
Lineages <- parallel::parLapply(cl,sample(basal_indexes),getLineage) %>% rbindlist
parallel::stopCluster(cl)
rm('cl')

ids <- Lineages[,sort(table(lineage_id),decreasing = T)]
i=0
i=i+1
Lineages[lineage_id==as.numeric(names(ids[i]))]



# ll <- getLineage(2693)
# Filtering lineages ------------------------------------------------------

### need to: 
#1) getLineageStat - F-stat from phylofactor-esque approximation
#2) filterIncompatibleLineages - remove lineages incompatible w/ F-maximizing ones
#3) fillLineageTable - get stat, P for NA's... will keep rank, rc-index NA'd
### May want to fill lineage table first - adding rc_index will put new edges on RCmap
### and allow fast ancestor/descendant computation for coarse-graining & getLineageStat.


#### We need to partition elements of X into groups given by edge-pair descendants.
#### This will give us a block structure across X and we can use within-block means
#### and run an ANOVA.

#### Ideally, it would be great to produce low-rank Xhat=wDv'
#### where D is a diagonal matrix of scores. However,
#### this requires orthonormal bases w, v to get D=w'Xv thus Xhat=w(w'Xv)v'
#### Our current bases - W,V - are not orthonormal. Is this really a requirement?


# getLineageStat ----------------------------------------------------------

### This requires efficient ways of finding the indexes of X which
### correspond to a given edge pair.

### We'll find a set of boxes in each lineage - starting with the entire lineage itself
### and that set of boxes can be refined. 

Xdt <- data.table('i'=rep(1:nrow(X),times=ncol(X)),
                     'j'=rep(1:ncol(X),each=nrow(X)),
                     'x'=c(X),
                  'lineage'=NA)

### lineages will be defined by row-tree tips


edgeTips <- function(tree){
  grps <- phylofactor::getPhyloGroups(tree)
  mat <- lapply(grps,FUN=function(gg)
                data.table('min'=min(gg[[1]]),'max'=max(gg[[1]])))
  mat <- rbindlist(mat)
  mat[,edge:=1:.N]
  setkey(mat,edge)
  mat[,tip:=(min==max)]
  return(mat)
}

colEdgeTips <- edgeTips(col.tree)
rowEdgeTips <- edgeTips(row.tree)

rowEdgeTips[edge==914,(max-min)] ##this edge carves out a very large lineage


lineageBoxes <- function(lineage,X,colEdgeTips,rowEdgeTips){
  #### if we start at the finest lineage w/ highest row.edge, we can construct bins iteratively
  #### 
  setkey(lineage,row.edge,col.edge)
  rowset <- rowEdgeTips[edge==lineage$row.edge[1],seq(min,max)]
  colset <- colEdgeTips[edge==lineage$col.edge[1],seq(min,max)]
  
  Xdt <- data.table('i'=rep(rowset,times=ncol(X)),
                    'j'=rep(1:ncol(X),each=length(rowset)),
                    'x'=c(X[rowset,]),
                    'lineage_id'=unique(lineage$lineage_id),
                    'box'=0)
  
  setkey(lineage,row.edge)
  setkey(Xdt,i,j)
  
  #### rowset will be continuously whittled down
  #### colset will need to be iterated.
  mbox=0
  for (nn in 1:nrow(lineage)){
    mbox=mbox+1
    if (nn==1){
      Xdt[j %in% colset,box:=mbox]
    } else {
      spp <- rowEdgeTips[edge==lineage$row.edge[nn],seq(min,max)]
      cols <- colEdgeTips[edge==lineage$col.edge[nn],seq(min,max)]
      bx <- Xdt[i %in% spp & box!=0,unique(box)] # The first box is the full set of species in the continents in which they're not found
      
      Xdt[i %in% spp & box!=0,box:=mbox]
      mbox=mbox+1
      Xdt[i %in% spp & j %in% cols,box:=mbox]
    }
  }
  return(Xdt)
}


id <- Lineages$lineage_id[1]
lineage <- Lineages[lineage_id==id]
Xdt <- lineageBoxes(lineage,X,colEdgeTips,rowEdgeTips)



######## THIS TAKES A WHILE - PARALLELIZE 
lineage_Fstats <- function(Lineages,X,colEdgeTips,rowEdgeTips){
  lineage_ids <- unique(Lineages$lineage_id)
  fits <- data.table('lineage_id'=lineage_ids,
                     'F_stat'=0)
  for (id in lineage_ids){
    Xdt <- lineageBoxes(Lineages[lineage_id==id],X,colEdgeTips,rowEdgeTips)
    setkey(Xdt,box)
    Xdt[,box:=factor(box)]
    a <- aov(x~box,data=Xdt)
    fits[lineage_id==id,F_stat:=summary(a)[[1]]['box','F value']]
  }
  return(fits)
}


lineage_ids <- Lineages[,list(n=.N),by=lineage_id][n>1]$lineage_id
fits <- lineage_Fstats(Lineages,X,colEdgeTips,rowEdgeTips)
fits <- fits[order(F_stat,decreasing = T)]
ix=0
ix=ix+1
Lineages[lineage_id==fits$lineage_id[ix],]
edg <- Lineages[lineage_id==fits$lineage_id[ix],row.edge[1]]
row.tree$tip.label[rowEdgeTips[edge==edg,seq(min,max)]]

#### However, this requires 
# getFstat <- function(X,lineage_table,W,V){
#   w <- W[,lineage_table$row.node-nrow(W),drop=F]
#   v <- V[,lineage_table$col.node-nrow(V),drop=F]
##   Xhat <- w %*% (t(w) %*% X %*% v) %*% t(v)
#   Xhat <- lineagefit(lineage_table,X)
#   ix <- which(Xhat!=0) #indexes for blocks being predicted by wDv'
#   n <- length(ix)
#   rss <- sum((X-mean(X)-Xhat)[ix]^2)
#   ess <- sum(c(Xhat[ix])^2)
#   dfm <- ncol(w)+1
#   df0 <- n-(ncol(w)+1)
#   Fstat <- (ess/dfm)/(rss/df0)
#   return(Fstat)
# }

### NOTE: for column tree, we need two bases: 

LineageONBasis <- function(lineage,W,V){
  row.edges <- sort(lineage$row.edges,decreasing = F)
  
  
}




# Iterating F-stats over P-value thresholds -------------------------------

source('R/edge_dendromap_fcns.R')
rc_table <- edge_rc_table(X,row.tree,col.tree,maxPval=0.01)

# edgeMap for quick descendant calculation -------------------------------------------------------------------
rowEdgeMap <- makeEdgeMap(row.tree)
colEdgeMap <- makeEdgeMap(col.tree)

colEdgeTips <- edgeTips(col.tree)
rowEdgeTips <- edgeTips(row.tree)

# RC_map ------------------------------------------------------------------
row.edges <- 1:Nedge(row.tree)
col.edges <- 1:Nedge(col.tree)
Row_Descendants <- rowEdgeMap[,list(list(setdiff(seq(edge,edge+n),edge))),by=edge]$V1
Col_Descendants <- colEdgeMap[,list(list(setdiff(seq(edge,edge+n),edge))),by=edge]$V1
RCmap <- edgeRCMap(rc_table,Row_Descendants,Col_Descendants)


# initialize cluster ------------------------------------------------------
# 
# cl <- parallel::makeCluster(7)
# 
# parallel::clusterExport(cl,varlist=c('incompatible_descendants','clean_sisters','greedyColTree','getLineage','edge_ancestors',
#                                      'RCmap','rc_table','Row_Descendants','Col_Descendants','row.tree','col.tree'),
#                         envir = environment())
# parallel::clusterEvalQ(cl,library(data.table))
# parallel::clusterEvalQ(cl,library(magrittr))



# iterate over Pvals -----------------------------------------
pvals <- sort(unique(rc_table$P),decreasing = F)

min_pval <-  cbind(rc_table[match(RCmap$ancestor,rc_index),P],
                       rc_table[match(RCmap$descendant,rc_index),P]) %>% apply(1,max) %>% min
RCmap$max_P <- cbind(rc_table[match(RCmap$ancestor,rc_index),P],
                     rc_table[match(RCmap$descendant,rc_index),P]) %>% apply(1,max)

min_ix=min(which(pvals==min_pval))
n_pvals=20
coarse_ix=round(seq(min_ix,length(pvals),length.out = n_pvals))
Fstats <- numeric(n_pvals)

Lineages <- NULL
# cl <- NULL
t_start=Sys.time()
for (i in coarse_ix){
  # i=i+1
    p_thresh=pvals[i]
    ix_thresh=rc_table[P<=p_thresh,rc_index]
    rct <- rc_table[rc_index %in% ix_thresh]
    rcm <- RCmap[max_P<=p_thresh]

    basal_indexes <- rcm[ancestor %in% ix_thresh & descendant %in% ix_thresh,setdiff(unique(ancestor),unique(descendant))]
    desc_count <- rcm[ancestor %in% ix_thresh & descendant %in% ix_thresh,list(n=.N),by=ancestor]
    basal_ixs_with_descendants <- intersect(basal_indexes,desc_count[n>1,ancestor])
    # basal_indexes
    # if (length(basal_ixs_with_descendants>0) & is.null(Lineages)){
      Lineages <- lapply(basal_ixs_with_descendants,getLineage,rc_table=rct,RCmap=rcm) %>% rbindlist # method=method
      Lineages[,lineage_size:=.N,by=lineage_id]
      Lineages <- Lineages[lineage_size>1]

      if (length(unique(Lineages$lineage_id))>1){ ### ensure lineages don't conflict, filter accordingly
        stats <- lineage_stats(Lineages,X,colEdgeTips,rowEdgeTips)
        stats <- stats[order(F_stat,decreasing = T)]
        stats <-  filter_stats(stats,Lineages,rct,rcm,Row_Descendants)
        Lineages <- Lineages[lineage_id %in% stats$lineage_id]
      }

      if (nrow(Lineages)>0){
        Fstats[match(i,coarse_ix)] <- global_Fstat(Lineages,X,colEdgeTips,rowEdgeTips)
      }
      
      
    # } else {
    #   new_indexes <- rct[P>pvals[match(i,coarse_ix)-1] & P<=p_thresh,rc_index]
    #   new_row.edges <- rct[rc_index %in% new_indexes,row.edge]
    #   # old_indexes <- setdiff(ix_thresh,new_indexes)
    #   ## easy way out: re-run all basal_ix in row.tree ancestor/descendant paths of new_ix edges
    #   row_ancs <- sapply(new_row.edges,edge_ancestors,tree=row.tree)
    #   row_descs <- rowEdgeMap[edge %in% new_row.edges,list(seq(edge,edge+n)),by=edge]$V1 %>% unique
    # 
    #   affected_ix <- rct[row.edge %in% c(row_ancs,row_descs,new_row.edges),rc_index]
    #   affected_basal_ix <- intersect(affected_ix,basal_ixs_with_descendants)
    # 
    #   if (length(affected_basal_ix)>0){
    # 
    #     if (!is.null(cl) & length(affected_basal_ix)>1){
    #       updated_lineages <- parallel::parLapply(cl,sample(affected_basal_ix),getLineage) %>% rbindlist
    #     } else {
    #       updated_lineages <- lapply(affected_basal_ix,getLineage,rc_table=rct,RCmap=rcm) %>% rbindlist
    #     }
    #     updated_lineages[,lineage_size:=.N,by=lineage_id]
    # 
    #     ## only recompute stats & filter affected lineages
    #     stats <- lineage_stats(updated_lineages,X,colEdgeTips,rowEdgeTips)
    #     stats <- stats[order(F_stat,decreasing = T)]
    #     stats <-  filter_stats(stats,updated_lineages,rct,rcm,Row_Descendants)
    #     updated_lineages <- updated_lineages[lineage_id %in% stats$lineage_id]
    # 
    # 
    # 
    #     affected_lineages <- Lineages[rc_index %in% affected_ix,unique(lineage_id)]
    # 
    #     Lineages <- rbind(Lineages[!lineage_id %in% affected_lineages],updated_lineages)
    #     Lineages[,lineage_size:=.N,by=lineage_id]
    #     Lineages <- Lineages[lineage_size>1]
    #     if (any(duplicated(Lineages[,c('row.edge','col.edge')]))){
    #       ## this can happen if e.g. for a given row edge, col.edge 8 is a descenant of another basal row.edge with col.edge6, but
    #       ## the same row.edge has significant col.edge 6 association and both descendants are col.edge 8, causing clean-up to assign
    #       ## col.edge 8 to the same row.edge, producing (i,8,rc_ix=1), (i,8,rc_ix=2)
    #       ### To resolve these disputes, we'll start off removing
    #       # stop('duplicated stuff!')
    # 
    #     }
    #     if (nrow(Lineages)>0){
    #       if (length(unique(Lineages$lineage_id))>1){
    #         Fstats[match(i,coarse_ix)] <- global_Fstat(Lineages,X,colEdgetips,rowEdgeTips)
    #       } else {
    #         Fstats[match(i,coarse_ix)] <- stats$F_stat
    #       }
    #     }
    #   } else { ## else: the new_index is a singleton edge and we don't modify the Lineages table
    #     Fstats[match(i,coarse_ix)] <- Fstats[match(i,coarse_ix)-1]
    #   }
    # 
    # }
}
t_stop=Sys.time()

t_stop-t_start
parallel::stopCluster(cl)
rm('cl')

continent_map <- data.table('col.edge'=1:8,
                            'continent'=c('Laurasia','Eurasia','NAmerica',
                                          'Gondwana','Australia','SA/Africa','Africa','SouthAmerica'))

save(list=ls(),file='data/edendromap_birds_workspace')

plot(pvals[coarse_ix],Fstats,log='x')

### next: we'll need to reduce the time and obtain the same result.

Lineages[row.edge==1] ## did we get rheas?

# stop cluster ------------------------------------------------------------


Lineages[,n:=.N,by=lineage_id]
Lineages <- Lineages[n>1]


# Filtering lineages ------------------------------------------------------

colEdgeTips <- edgeTips(col.tree)
rowEdgeTips <- edgeTips(row.tree)

### get fit stats
stats <- lineage_stats(Lineages,X,colEdgeTips,rowEdgeTips)
stats <- stats[order(F_stat,decreasing = T)]
stats <-  filter_stats(stats,Lineages,rc_table,RCmap,Row_Descendants)

Lineages <- Lineages[lineage_id %in% stats$lineage_id]