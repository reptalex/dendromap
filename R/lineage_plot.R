#' Focused of particular lineage
#' @export
#' @param dm dendromap object
#' @param lineage lineage number in \code{dm$Lineage} for plotting
#' @param id Alternative to \code{lineage}, \code{dm$lineage_id} for plotting
#' @param ... additional input arguments for \code{\link{plot.dendromap}}
#' @examples
lineage_plot <- function(dm,lineage=NULL,id=NULL,...){
  ### trim row.tree to only lineage; rename/follow edges for right labelling
  ### plot subset of data
  
  if (is.null(id)){
    L <- dm$Lineages[Lineage==lineage]
  } else {
    L <- dm$Lineages[lineage_id==id]
  }
  
  
  basal_edge <- L[,min(row.edge)]
  desc_edges <- dm$rowDescendants[[basal_edge]]
  tips <- phangorn::Descendants(dm$row.tree,dm$row.tree$edge[basal_edge,2],'tips')[[1]]
  rt <- ape::drop.tip(dm$row.tree,setdiff(dm$row.tree$tip.label,dm$row.tree$tip.label[tips]))
  rownames(rt$edge) <- desc_edges
  
  L[row.edge!=basal_edge,row.node:=rt$edge[as.character(row.edge),2]]
  L[row.edge==basal_edge,row.node:=(length(rt$tip.label)+1)]
  
  temp_dm <- list('Data'=dm$Data[row.tree$tip.label[tips],],
                  'Lineages'=L,
                  'col.tree'=dm$col.tree,
                  'row.tree'=rt)
  class(temp_dm) <- 'dendromap'
  pl <- plot.dendromap(temp_dm,...)
  return(pl)
}