#' Make table of {row tree,column tree} edge pair statistics
#' @export
#' @param X dataset whose rows correspond exactly to \code{row.tree$tip.label} and columns correspond to \code{col.tree$tip.label}. NA values not allowed.
#' @param row.tree \code{phylo} class object
#' @param col.tree \code{phylo} class object
#' @param maxPval Cutoff - remove all P-values greater than \code{maxPval}
#' @examples
#' ## none yet
make_rc_table <- function(X,row.tree,col.tree,maxPval=0.01){
  m <- nrow(X)
  n <- ncol(X)
  W <- edge_contrasts(row.tree)
  V <- edge_contrasts(col.tree)
  
  U <- t(W) %*% X %*% V
  Unull <- t(W) %*% X[sample(m),sample(n)] %*% V
  cdf <- ecdf(c(Unull))
  
  
  rc_table=which(abs(U)>0) %>% arrayInd(.dim = c(nrow(U),ncol(U))) %>% as.data.table
  names(rc_table) <- c('row.edge','col.edge')
  rc_table$stat <- U[abs(U)>0]
  setkey(rc_table,row.edge,col.edge,stat)
  
  rc_table[,rank:=rank(-abs(stat))]
  rc_table[,P:=1-cdf(stat)]
  rc_table <- rc_table[stat>0]
  setkey(rc_table,stat)
  
  y <- log(c(Unull[Unull>0]))
  if (any(rc_table$P==0)){
    nn=sum(rc_table$P==0)
    base::cat(paste('\n',nn,' P-values were 0. Will estimate tail probabilities assuming log(stat^2)~rnorm(mu,sd)',sep=''))
    y <- y[y>-20]
    mu <- mean(y)
    sig <- sd(y)
    min.P <- rc_table[P>0,min(P)]
    stat.min <- max(log(rc_table[P>0][P==min(P),stat]))
    estimate_tail <- function(y1,mu,sig,ymin,pmin) (1-pnorm(y1,mu,sig))/(1-pnorm(ymin,mu,sig))*pmin
    rc_table[P==0,P:=estimate_tail(log(stat),mu,sig,stat.min,min.P)]
  }
  
  rc_table <- rc_table[P<=maxPval]
  rc_table[,rc_index:=1:.N]
  
  return(rc_table)
}