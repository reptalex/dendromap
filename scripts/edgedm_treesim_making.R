rm(list=ls())
load('data/birds/bird_dendromap_workspace')
source('R/edge_dendromap_fcns.R')
library(phylofactor)
library(dendromap)

set.seed(1)
row.tree <- rtree(50)

par(mfrow=c(1,2))
plot(row.tree)
edgelabels()
plot(col.tree)
edgelabels()

# edgeDepths <- function(tree){
#   x <- data.table('edge'=1:nrow(tree$edge))
# }

treeSim <- function(n,row.tree,col.tree,row.edge=NULL,col.edge=NULL,row.edge=NULL
                    fix.col.edge=FALSE,depth_matching=FALSE){
  
  if (is.null(col.edge)){
    fix.col.edge <- FALSE
  }
  if (!is.null(row.edge)){
    n=1
  }
  
  if (depth_matching){ ## unfinished
  }
  
  row.edges <- 1:nrow(row.tree$edge)
  interior_row_edges <- which(row.tree$edge[,2]>length(row.tree$tip.label))
  col.edges <- 1:nrow(col.tree$edge)
  interior_col_edges <- which(row.tree$edge[,2]>length(col.tree$tip.label))
  
  rowEdgeMap <- makeEdgeMap(row.tree)
  colEdgeMap <- makeEdgeMap(col.tree)
  Row_Descendants <- rowEdgeMap[,list(list(setdiff(seq(edge,edge+n),edge))),by=edge]$V1
  Col_Descendants <- colEdgeMap[,list(list(setdiff(seq(edge,edge+n),edge))),by=edge]$V1
  
  for (i in 1:n){
    
    row.edge <- sample(interior_row_edges,1)
    col.edge <- sample(interior_col_edges,1)
    
    lineage <- data.table('row.edge'=row.edge,'col.edge'=col.edge)
    
    
    
    done=F
    while(!done){
      
    }
    
    
    
  }
  
}