#' Random dendromap object
#' @export
#' @param m number of tips in \code{row.tree}
#' @param n number of tips in \code{col.tree}
#' @param row.tree optional phylo class input, instead of m.
#' @param col.tree optional phylo class input, instead of n.
#' @param lineage_sizes either single integer or integer vector of length \code{n_lineages}. Lineages will be less than or equal to \code{lineage_sizes}
#' @param lambda positive numeric rate parameter for \code{\link{rpois}} generating lineage sizes
#' @param n_lineages integer number of n_lineages to compute
#' @param min.lineage.size minimum lineage size - 
rdendromap <- function(m=NULL,n=NULL,row.tree=NULL,col.tree=NULL,lineage_sizes=NULL,
                       lambda=3,n_lineages=1,min.lineage.size=2,simulate.dataset=FALSE,...){
  
  if (is.null(m) & is.null(row.tree)){
    stop('Must input either m or row.tree')
  }
  if (is.null(n) & is.null(col.tree)){
    stop('Must input either n or col.tree')
  }

  row.tree <- NULL
  col.tree <- NULL
  if (is.null(row.tree)){
    row.tree <- ape::rtree(m)
  }
  if (is.null(col.tree)){
    col.tree <- ape::rtree(n)
  }
  
  # row_height <- ape::node.height(row.tree)
  # col_height <- ape::node.height(col.tree)
  
  rowEdgeMap <- edge_registry(row.tree)
  colEdgeMap <- edge_registry(col.tree)
  rowEdgeTips <- edge_tips(row.tree)
  colEdgeTips <- edge_tips(col.tree)
  rowDescendants <- rowEdgeMap[,list(list(setdiff(seq(edge,edge+n),edge))),by=edge]$V1
  colDescendants <- colEdgeMap[,list(list(setdiff(seq(edge,edge+n),edge))),by=edge]$V1
  if (is.null(lineage_sizes)){
    lineage_sizes <- rep(0,n_lineages)
    small <- lineage_sizes<min.lineage.size
    while (any(small)){
      lineage_sizes[small] <- rpois(sum(small),lambda)
      small <- lineage_sizes<min.lineage.size
    }
  } else {
    if (length(lineage_sizes)==1){
      lineage_sizes <- rep(lineage_sizes,n_lineages)
    } else {
      if (length(lineage_sizes)!=n_lineages){
        stop('Input lineage_sizes must either be integer or vector of length equal to n_lineages')
      }
    }
  }
  remaining_row_edges <- which(rowEdgeMap[,n>=(min.lineage.size-1)])
  Lineages <- NULL
  for (nn in 1:n_lineages){
    basal_row_edge <- sample(rowEdgeMap[edge %in% remaining_row_edges & n>=lineage_sizes[nn],edge],size=1)
    basal_col_edge <- sample(colEdgeMap[n>1,edge],size=1)
    
    lineage <- grow_lineage(basal_row_edge,basal_col_edge,
                            row.tree,col.tree,
                            rowEdgeMap,colEdgeMap,lineage_sizes[nn],
                            rowDescendants,colDescendants)
    lineage[,lineage_id:=nn]
    Lineages <- rbind(Lineages,lineage)
    
    remaining_row_edges <- setdiff(remaining_row_edges,rowDescendants[[basal_row_edge]])
    remaining_row_edges <- setdiff(remaining_row_edges,edge_ancestors(basal_row_edge,row.tree))
    if (length(remaining_row_edges)==0){
      break
    }
  }
  
  edge2node <- function(edges,tree) tree$edge[edges,2]
  Lineages[,row.node:=edge2node(row.edge,row.tree)]
  Lineages[,col.node:=edge2node(col.edge,col.tree)]
  
  object <- list('Lineages'=Lineages,
                 'Data'=NA,
                 'row.tree'=row.tree,
                 'col.tree'=col.tree,
                 'colEdgeTips'=colEdgeTips,
                 'rowEdgeTips'=rowEdgeTips,
                 'rowDescendants'=rowDescendants,
                 'colDescendants'=colDescendants)
  
  class(object) <- 'dendromap'
  if (simulate.dataset){
    object$Data <- simulate(object)
  }
  return(object)
}
