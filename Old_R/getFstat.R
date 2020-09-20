#' Calculate F-statistic for low-rank approximation X=WDV'+E for lineages in lineage_table
#' @export
#' @param X data matrix
#' @param lineage_table output from e.g. \code{\link{findLineageTable}}
#' @param W row tree-basis, \code{\link{treeBasis}}
#' @param V col tree-basis, \code{\link{treeBasis}}
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
