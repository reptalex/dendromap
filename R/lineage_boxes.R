#' Coarse-grain dataset, X, based on discovered lineage
#' @export
#' @param lineage see \code{\link{get_lineage}}
#' @param X dataset input to \code{\link{dendromap}} from which lineages were made
#' @param colEdgeTips \code{\link{edge_tips}} for col.tree
#' @param rowEdgeTips \code{\link{edge_tips}} for row.tree
#' @examples
lineage_boxes <- function(lineage,X,colEdgeTips,rowEdgeTips){
  #### if we start at the finest lineage w/ highest row.edge, we can construct bins iteratively
  #### 
  setkey(lineage,row.edge,col.edge)
  rowset <- rowEdgeTips[edge==lineage$row.edge[1],seq(min,max)]
  colset <- colEdgeTips[edge==lineage$col.edge[1],seq(min,max)]
  
  Xdt <- data.table('i'=rep(rowset,times=ncol(X)),
                    'j'=rep(1:ncol(X),each=length(rowset)),
                    'x'=c(X[rowset,]),
                    'lineage_id'=unique(lineage$lineage_id),
                    'box'=0)
  
  setkey(lineage,row.edge)
  setkey(Xdt,i,j)
  
  #### rowset will be continuously whittled down
  #### colset will need to be iterated.
  mbox=0
  for (nn in 1:nrow(lineage)){
    mbox=mbox+1
    if (nn==1){
      Xdt[j %in% colset,box:=mbox]
    } else {
      spp <- rowEdgeTips[edge==lineage$row.edge[nn],seq(min,max)]
      cols <- colEdgeTips[edge==lineage$col.edge[nn],seq(min,max)]
      bx <- Xdt[i %in% spp & box!=0,unique(box)] # The first box is the full set of species in the continents in which they're not found
      
      Xdt[i %in% spp & box!=0,box:=mbox]
      mbox=mbox+1
      Xdt[i %in% spp & j %in% cols,box:=mbox]
    }
  }
  return(Xdt)
}