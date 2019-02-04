#' Map events on column tree to nodes on row tree
#' @export
#' @param row.tree \code{phylo} class object
#' @param col.tree \code{phylo} class object
#' @param row.node integer. optional node for starting treeMap. Must be a node in \code{row.tree}
#' @param col.node integer. optional node for starting treeMap. Must be a node in \code{col.tree}
#' @param prob.row probability, between 0 and 1, that a descendant node is manifested in the tree map.
#' @param prob.col same as \code{prob.row} but for nodes in column tree.
#' @examples
#' set.seed(1)
#' library(ape)
#' row.tree <- rtree(100)
#' col.tree <- rtree(10)
#' S <- treeMap(row.tree,col.tree,101,11,0.2,1)
#' pathViz(S,col.tree,row.tree)

treeMap <- function(row.tree,col.tree,row.node=NULL,col.node=NULL,
                    prob.row=1,prob.col=1){
  if (is.null(row.node)){
    row.node <- ape::Ntip(row.tree)+1
  }
  if (is.null(col.node)){
    col.node <- ape::Ntip(col.tree)+1
  }
  
  Path <- data.table('row.node'=row.node,
                     'col.node'=col.node,
                     'orientation'=sign(rnorm(1)),
                     'terminated'=F)
  
  row.nb <- nodeBank(ape::Ntip(row.tree)+1,row.tree,prob.row)
  col.nb <- nodeBank(ape::Ntip(col.tree)+1,col.tree,prob.col)
  while (any(!Path$terminated)){
    ix <- min(which(!Path$terminated))
    ch.row <- labelledChildren(Path$row.node[ix],row.tree)
    ch.col <- labelledChildren(Path$col.node[ix],col.tree)
    orientation <- Path$orientation[ix]
    if (length(ch.col)>0 & length(ch.row)>0){
      colA <- manifestChildren(ch.col['a'],col.tree,col.nb)
      rowA <- manifestChildren(ch.row['a'],row.tree,row.nb)
      colB <- manifestChildren(ch.col['b'],col.tree,col.nb)
      rowB <- manifestChildren(ch.row['b'],row.tree,row.nb)
      Aset <- list('a'=colA,'b'=colB)
      
      if (orientation>0){
        A <- (length(rowA)>0 & length(colA)>0)
        B <- (length(rowB)>0 & length(colB)>0)
      } else {
        A <- (length(rowA)>0 & length(colB)>0)
        B <- (length(rowB)>0 & length(colA)>0)
      }
      if (A | B){
        if (A){
          rA <- data.table('row.node'=rowA,
                           'side'='a',
                           'x'=runif(length(rowA)))
          side <- c('a','b')[as.numeric(orientation<0)+1]
          cA <- data.table('col.node'=Aset[[side]],
                           'side'=side,
                           'x'=runif(length(Aset[[side]])))
        } 
        
        if (B){
          rB <- data.table('row.node'=rowB,
                           'side'='b',
                           'x'=runif(length(rowB)))
          side <- c('a','b')[as.numeric(orientation>0)+1]
          cB <- data.table('col.node'=Aset[[side]],
                           'side'=side,
                           'x'=runif(length(Aset[[side]])))
        }
        
        if (A & B){
          row.children <- rbind(rA,rB)
          col.children <- rbind(cA,cB)
        } else {
          if (A & !B){
            row.children <- rA
            col.children <- cA
          } else {
            row.children <- rB
            col.children <- cB
          }
        }
        if (orientation<0){
          sw <- c('a','b')
          names(sw) <- c('b','a')
          col.children[,side:=sw[side]]
        }
        
        row.children[,x:=as.numeric(order(x,decreasing = T)),by=side]
        col.children[,x:=as.numeric(order(x,decreasing = T)),by=side]
        setkey(row.children,side,x)
        setkey(col.children,side,x)
        
        map <- row.children[col.children,nomatch=0]
        map[,orientation:=sample(c(1,-1),.N,T,c(.5,.5))]
        map[,terminated:= (row.node<=ape::Ntip(row.tree) | col.node<=ape::Ntip(col.tree))]
        map[row.node %in% row.nb[terminal==TRUE,node] |
              col.node %in% col.nb[terminal==TRUE,node]]$terminated <- TRUE
        Path <- rbind(Path,map[,c('row.node','col.node','orientation','terminated')])
        Path$terminated[ix] <- TRUE
      } else {
        Path$terminated[ix] <- TRUE
      }
    } else {
      Path$terminated[ix] <- TRUE
    }
  }
  return(Path)
}
