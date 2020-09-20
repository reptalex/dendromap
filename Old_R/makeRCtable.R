#' make table of significant row-column node pairs
#' @export
#' @param N matrix whose row names are in \code{row.tree$tip.label} and column names are in \code{col.tree$tip.label}
#' @param row.tree phylo class object
#' @param col.tree phylo class object
#' @param W \code(treeBasis(row.tree))
#' @param V \code(treeBasis(col.tree))
#' @param n_sim number of simulatd null datasets (used to compuate approximate P values for r-c node pairs)
#' @param cl cluster with \code{library(dendromap)} loaded
#' @param discard_contingency_zeros logical. If TRUE, will discard all rc pairs for which either (a) a sister clade is always zero across sites or (b) both clades are zero within a site.
makeRCtable <- function(N,row.tree,col.tree,W=NULL,V=NULL,n_sim=NULL,cl=NULL,discard_contingency_zeros=FALSE){
  
  if (is.null(W)){
    W <- treeBasis(row.tree)
  }
  if (is.null(V)){
    V <- treeBasis(col.tree)
  }
  U <- t(W) %*% as.matrix(N) %*% V
  row.nodes <- sapply(rownames(U),strsplit,'_') %>% sapply(getElement,2) %>% as.numeric
  col.nodes <- sapply(colnames(U),strsplit,'_') %>% sapply(getElement,2) %>% as.numeric
  
  if (discard_contingency_zeros){
    ix <- idComparables(N,W,V,cl)
  } else {
    ix <- data.table('i'=rep(1:ncol(W),times=m),
                     'j'=rep(1:ncol(V),each=n),
                     'comparable'=TRUE)
  }
  
  U[as.matrix(ix[comparable==FALSE,c('i','j')])] <- NA
  
  rc_table <- data.table(expand.grid('row.node'=row.nodes,
                                     'col.node'=col.nodes),
                         'stat'=c(U),key=c('row.node','col.node'))
  rc_table <- rc_table[!is.na(stat)]
  gc()
  if (is.null(n_sim)){
    n_sim <- ceiling(10000/(nrow(rc_table)))
  }
  
  rsim <- function(x,N,W,V,cl,discard_contingency_zeros){
    Nnull <- as.matrix(N[sample(nrow(N)),sample(ncol(N))])
    if (discard_contingency_zeros){
      base::cat('\nIdentifying comparable species|species x continent|continent splits in null dataset')
      ix <- idComparables(Nnull,W,V,cl,verbose=FALSE)
    } else {
      ix <- data.table('i'=rep(1:ncol(W),times=m),
                       'j'=rep(1:ncol(V),each=n),
                       'comparable'=TRUE)
    }
    Unull <- t(W) %*% Nnull %*% V
    Unull[as.matrix(ix[comparable==FALSE,c('i','j')])] <- NA
    log((c(Unull[!is.na(Unull)]))^2) %>%
      return()
  }
  
  y <- sapply(1:n_sim,rsim,N,W,V,cl,discard_contingency_zeros) %>% unlist
  null_cdf <-  ecdf(y[!is.infinite(y)])
  rc_table[,P:=1-null_cdf(log(stat^2))]
  if (any(rc_table$P==0)){
    nn=sum(rc_table$P==0)
    base::cat(paste('\n',nn,' P-values were 0. Will estimate tail probabilities assuming log(stat^2)~rnorm(mu,sd)',sep=''))
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
