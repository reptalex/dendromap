#' find best sequence of node-pairs in the correct orientation as the column tree
#' @param RowCol vector of length two containing indexes \code{(i,j)} mapping the row tree node represented in row \code{i} of \code{U} and column tree node represented in column \code{j} of \code{U}
#' @param U matrix whose entry in row i column j is the bilinear map from row-tree's i'th node and the column tree's j'th node
#' @param row.tree phylo class object
#' @param col.tree phylo class object
#' @param row.nodemap nodemap produced in \code{\link{makeNodeMap}} 
#' @param col.nodemap nodemap produced in \code{\link{makeNodeMap}}
findNodeSeq <- function(RowCol,U.=U,row.tree.=row.tree,col.tree.=col.ree,row.nodemap.=row.nodemap,col.nodemap.=col.nodemap){
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
