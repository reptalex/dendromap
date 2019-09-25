
#' dendromap - find cophylogenetic patterns in a dataset
#' @export
#' @param X matrix whose rownames are in \code{row.tree$tip.label} and colnames are in \code{col.tree$tip.label}
#' @param row.tree phylo class object. Polytomies will be ignored.
#' @param col.tree phylo class object. Polytomies will be ignored.
#' @param fdr_threshold number between 0 or 1, the upper bound approximate fdr of row-col node-pair statistics for inclusion in our search
#' @param Pval_threshold Threshold "significance" of a row-column node pair for consideration in lineages. If input, will override \code{fdr_threshold}
#' @param W optional \code{treeBasis(row.tree)} - must have colnames in the format of e.g. "node_51" for node 51.
#' @param V optional \code{treeBasis(col.tree)} - must have colnames in the format of e.g. "node_51" for node 51.
#' @param n_sim optional number of null datasets to simulate (via row & column shuffling) in order to generate approximate P-values in \code{\link{makeRCtable}}
#' @examples
#' library(dendromap)
#' set.seed(1)
#' m=1000
#' n=10
#' row.tree <- rtree(m)
#' col.tree <- rtree(n)
#' S <- treeSim(5,row.tree,col.tree,col.node=n+1)
#' eta <- S$W %*% (10*S$D) %*% t(S$V)
#' X <- eta+matrix(rnorm(m*n),nrow=m)
#' clrinv <- function(x) exp(x)/sum(exp(x))
#' rmlt <- function(p,lambda=5e3) rmultinom(1,rpois(1,lambda),prob = p)
#' N <- apply(X,2,clrinv) %>% apply(2,rmlt)
#' rownames(N) <- row.tree$tip.label
#' colnames(N) <- col.tree$tip.label
#' dm <- dendromap(N,row.tree,col.tree,W=S$W,V=S$V,Pval_threshold=0.1) 
#' 
#' plot.dendromap(S)
#' plot.dendromap(dm)
#' S$Lineages[,rc:=paste(row.node,col.node,sep='_')]
#' dm$Lineages[,rc:=paste(row.node,col.node,sep='_')]
#' sum(S$Lineages$rc %in% dm$Lineages$rc)/nrow(S$Lineages)      ### 78.6% of the real rc's were ID'd
#' sum(dm$Lineages$rc %in% S$Lineages$rc)/nrow(dm$Lineages)  ### at a 100% true positive rate

dendromap <- function(X,row.tree,col.tree,fdr_threshold=0.2,Pval_threshold=NULL,W=NULL,V=NULL,n_sim=NULL){
  
  base::cat(paste('Checking Data and tree compatibility'))
  ### Align dataset to trees
  if (!all(rownames(X) %in% row.tree$tip.label)){
    stop('There are rownames(X) not in row.tree$tip.label')
  } else {
    if (!all(row.tree$tip.label %in% rownames(X))){
      row.tree <- ape::drop.tip(row.tree,setdiff(rownames(X),row.tree$tip.label))
    }
    X <- X[row.tree$tip.label,]
  }
  if (!all(colnames(X) %in% col.tree$tip.label)){
    stop('There are colnames(X) not in col.tree$tip.label')
  } else {
    if (!all(col.tree$tip.label %in% colnames(X))){
      col.tree <- ape::drop.tip(col.tree,setdiff(colnames(X),col.tree$tip.label))
    }
    X <- X[,col.tree$tip.label]
  }

  base::cat(paste('\nMaking Nodemaps'))
  row.nodemap <- dendromap:::makeNodeMap(row.tree)
  col.nodemap <- dendromap:::makeNodeMap(col.tree)

  base::cat(paste('\nMaking RC table with',n_sim,'null simulations'))
  rc_table <- makeRCtable(X,row.tree,col.tree,W,V,n_sim)
  rc_table <- rc_table[P<max(P)]
  rc_table[,fdr:=p.adjust(P,'fdr')]
  if (!is.null(Pval_threshold)){
    rc_table <- rc_table[P<=Pval_threshold]
    base::cat(paste('\n',nrow(rc_table),' RCs had P<=Pval_threshold at Pval_threshold=',Pval_threshold,sep=''))
  } else {

    rc_table <- rc_table[fdr<=fdr_threshold]
    base::cat(paste('\n',nrow(rc_table),' RCs had fdr<fdr_threshold at fdr_threshold=',Pval_threshold,sep=''))
  }
  
  RCmap <- makeRCMap(rc_table,row.nodemap,col.nodemap)
  
  
  base::cat(paste('\nRCmap has',nrow(RCmap),'rows'))
  Lineages <- find_lineages(RCmap,rc_table,row.nodemap,col.nodemap)
  compute_score <- function(lineage,rc_table.=rc_table) rc_table[rc_index %in% lineage,-sum(log(P))]
  
  base::cat(paste('\nRC-RC Tree traversal found',length(Lineages),'sequences of RCs. \n If this number is large, joining sequences by finding cliques will take a long time.'))
  base::cat(paste('\nFiltering RC sequences into lineages'))
  
  i=0
  output <- NULL
  while (length(Lineages)>0){
    i=i+1

    scores <- sapply(Lineages,compute_score)
    winner <- which.max(scores)
    
    output_table <- rc_table[rc_index %in% Lineages[[winner]]]
    output_table[,Lineage:=i]
    output <- rbind(output,output_table)
    Lineages <- filter_winner(winner,Lineages,rc_table,row.nodemap)
  }
  output <- list('Lineages'=output,'Data'=X,'row.tree'=row.tree,'col.tree'=col.tree)
  class(output) <- 'dendromap'
  return(output)
}
