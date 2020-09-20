#' Predict dataset from dendromap object
#' @export
#' @param object \code{\link{dendromap}} object
#' @param ... additional arguments, currently none are used
#' @examples
predict.dendromap <- function(object,...){
  
  lineages <- split(object$Lineages,f = factor(object$Lineages$lineage_id)) %>% 
    lapply(lineage_boxes,object$Data,object$colEdgeTips,object$rowEdgeTips) %>% rbindlist
  
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
  Xdt <- matrix(Xdt$yhat,nrow=nrow(X),ncol=ncol(X),byrow=T)
  rownames(Xdt) <- rownames(object$Data)
  colnames(Xdt) <- colnames(object$Data)
  return(Xdt)
  
}