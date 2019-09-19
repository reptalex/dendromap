#' make table of significant row-column node pairs
#' @export
#' @param N matrix whose row names are in \code{row.tree$tip.label} and column names are in \code{col.tree$tip.label}
#' @param row.tree phylo class object
#' @param col.tree phylo class object
#' @param W \code(treeBasis(row.tree))
#' @param V \code(treeBasis(col.tree))
#' @param n_sim number of simulatd null datasets (used to compuate approximate P values for r-c node pairs)
makeRCtable <- function(N,row.tree,col.tree,W=NULL,V=NULL,n_sim=NULL){
  
  if (is.null(W)){
    W <- treeBasis(row.tree)
  }
  if (is.null(V)){
    V <- treeBasis(col.tree)
  }
  U <- t(W) %*% as.matrix(N) %*% V
  row.nodes <- sapply(rownames(U),strsplit,'_') %>% sapply(getElement,2) %>% as.numeric
  col.nodes <- sapply(colnames(U),strsplit,'_') %>% sapply(getElement,2) %>% as.numeric
  
  rc_table <- data.table(expand.grid('row.node'=row.nodes,
                                     'col.node'=col.nodes),
                         'stat'=c(U))
  gc()
  if (is.null(n_sim)){
    n_sim <- ceiling(10000/(nrow(N)*ncol(N)))
  }
  rsim <- function(x,N,W,V) abs(c(t(W) %*% as.matrix(N[sample(nrow(N)),sample(ncol(N))]) %*% V))
  null_cdf <- sapply(1:n_sim,rsim,N,W,V) %>% ecdf
  rc_table[,P:=1-null_cdf(abs(stat))]
  if (any(rc_table$P==0)){
    warning('Some P-values are 0. Will replace with 1/(n_sim*nrow(N)*ncol(N))')
    rc_table[P==0,P:=1/(2*n_sim*nrow(N)*ncol(N))]
  }
  rc_table[,rc_index:=1:.N]
  setkey(rc_table,row.node,col.node)
  rm(list=c('U','W','V'))
  gc()
  return(rc_table)
}
