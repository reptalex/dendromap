#' Extract global F statistic for several Lineages
#' @export
#' @param Lineages see \code{\link{get_lineages_and_stats}}
#' @param X dataset input to \code{\link{dendromap}} from which lineages were made
#' @param colEdgeTips \code{\link{edge_tips}} for col.tree
#' @param rowEdgeTips \code{\link{edge_tips}} for row.tree
#' @examples
global_Fstat <- function(Lineages,X,colEdgeTips,rowEdgeTips){
  lineages <- split(Lineages,f = factor(Lineages$lineage_id)) %>% 
    lapply(lineage_boxes,X,colEdgeTips,rowEdgeTips) %>% rbindlist
  
  Xdt <- data.table('i'=rep(1:nrow(X),times=ncol(X)),
                    'j'=rep(1:ncol(X),each=nrow(X)),
                    'x'=c(X))
  setkey(Xdt,i,j)
  setkey(lineages,i,j)
  Xdt <- lineages[,c('i','j','lineage_id','box')][Xdt]
  Xdt[is.na(lineage_id),lineage_id:=-1]
  Xdt[is.na(box),box:=-1]
  Xdt[,group:=paste(lineage_id,box,sep='_')]
  Xdt[,yhat:=mean(x,na.rm=T),by=group]
  a=aov(x~box,data=Xdt) %>% summary
  return(a[[1]]['box','F value'])
}