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
                         'stat'=c(U),key=c('row.node','col.node'))
  gc()
  if (is.null(n_sim)){
    n_sim <- ceiling(10000/(nrow(N)*ncol(N)))
  }
  rsim <- function(x,N,W,V) log((c(t(W) %*% as.matrix(N[sample(nrow(N)),sample(ncol(N))]) %*% V))^2)
  y <- sapply(1:n_sim,rsim,N,W,V)
  null_cdf <-  ecdf(y)
  rc_table[,P:=1-null_cdf(log(stat^2))]
  if (any(rc_table$P==0)){
    warning('Some P-values were 0. Will estimate tail probabilities assuming log(stat^2)~rnorm(mu,sd). You may want to consider manually increasing n_sim for more accurate null distribution')
    y <- y[y>-20]
    mu <- mean(y)
    sig <- sd(y)
    min.P <- rc_table[P>0,min(P)]
    stat.min <- max(log(rc_table[P>0][P==min(P),stat]^2))
    estimate_tail <- function(y1,mu,sig,ymin,pmin) (1-pnorm(y1,mu,sig))/(1-pnorm(ymin,mu,sig))*pmin
    rc_table[P==0,P:=estimate_tail(log(stat^2),mu,sig,stat.min,min.P)]
  }
  rc_table[,rc_index:=1:nrow(rc_table)]
  rm(list=c('U','W','V'))
  gc()
  return(rc_table)
}
