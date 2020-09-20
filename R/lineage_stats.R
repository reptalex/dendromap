#' Extract F statistics for all lineages
#' @export
#' @param Lineages see \code{\link{get_lineage}}
#' @param X dataset input to \code{\link{dendromap}} from which lineages were made
#' @param colEdgeTips \code{\link{edge_tips}} for col.tree
#' @param rowEdgeTips \code{\link{edge_tips}} for row.tree
#' @param cl cluster initialized internally with \code{\link{dendromap}} 
#' @examples
lineage_stats <- function(Lineages,X,colEdgeTips,rowEdgeTips,cl=NULL){
  lineage_ids <- unique(Lineages$lineage_id)
  fits <- data.table('lineage_id'=lineage_ids,
                     'F_stat'=0)
  if (is.null(cl)){
    for (id in lineage_ids){
      Xdt <- lineage_boxes(Lineages[lineage_id==id],X,colEdgeTips,rowEdgeTips)
      setkey(Xdt,box)
      Xdt[,box:=factor(box)]
      a <- aov(x~box,data=Xdt)
      ss <-summary(a)[[1]]['box',] %>% unlist
      fits[lineage_id==id,F_stat:=ss['F value']]
      fits[lineage_id==id,msq:=ss['Mean Sq']]
      fits[lineage_id==id,ssq:=ss['Sum Sq']]
    }
  } else {
    fits <- parallel::parLapply(cl,lineage_ids,lstat,Lineages=Lineages) %>% rbindlist
  }
  return(fits)
}