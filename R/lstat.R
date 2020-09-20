#' Extract F statistic for a given lineage
#' @export
#' @param id value found in \code{Lineages$lineage_id}
#' @param Lineages see \code{\link{get_lineage}}
#' @param X dataset input to \code{\link{dendromap}} from which lineages were made
#' @param colEdgeTips \code{\link{edge_tips}} for col.tree
#' @param rowEdgeTips \code{\link{edge_tips}} for row.tree
#' @examples
lstat <- function(id,Lineages,X.=X,colEdgeTips.=colEdgeTips,rowEdgeTips.=rowEdgeTips){
  Xdt <- lineage_boxes(Lineages[lineage_id==id],X,colEdgeTips,rowEdgeTips)
  setkey(Xdt,box)
  Xdt[,box:=factor(box)]
  a <- aov(x~box,data=Xdt)
  ss <-summary(a)[[1]]['box',] %>% unlist
  output <- data.table('lineage_id'=id,'F_stat'=ss['F value'],'msq'=ss['Mean Sq'],'ssq'=ss['Sum Sq'])
  return(output)
}