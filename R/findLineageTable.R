#' Find table of lineages given rc_table subset by a P-value threshold
#' @export
#' @param rc_tbl table from \code{\link{makeRCtable}}
#' @param row.nodemap row tree nodemap from \code{\link{makeNodeMap}}
#' @param col.nodemap col tree nodemap from \code{\link{makeNodeMap}}
#' @param cl optional input cluster with package \code{dendromap} functions
findLineageTable <- function(rc_tbl,row.nodemap,col.nodemap,cl){
  if (nrow(rc_tbl)>1){
    row.nodes <- unique(rc_tbl$row.node)
    col.nodes <- unique(rc_tbl$col.node)
    Row_Descendants <- lapply(row.nodes,getIndexSets,row.nodemap) %>%
      lapply(FUN=function(x,a) lapply(x,intersect,a),a=row.nodes)
    Col_Descendants <- lapply(col.nodes,getIndexSets,col.nodemap) %>%
      lapply(FUN=function(x,a) lapply(x,intersect,a),a=col.nodes)
    
    names(Row_Descendants) <- row.nodes
    names(Col_Descendants) <- col.nodes
    RCmap <- tryCatch(makeRCMap(rc_tbl,Row_Descendants,Col_Descendants),error=function(e) NULL)
    if (is.null(RCmap) | all(RCmap$terminal)){
      rc_tbl[,Lineage:=1:.N]
      return(rc_tbl)
    } else {
      # base::cat(paste('\nRCmap has',nrow(RCmap),'rows'))
      Lineages <- find_lineages(RCmap,rc_tbl,Row_Descendants,Col_Descendants,cl)
      compute_score <- function(lineage,rc_tbl.=rc_tbl) rc_tbl[rc_index %in% lineage,-sum(log(P))]
      
      # base::cat(paste('\nRC-RC Tree traversal found',length(Lineages),'sequences of RCs. \n If this number is large, joining sequences by finding cliques will take a long time.'))
      # base::cat(paste('\nFiltering RC sequences into lineages'))
      
      i=0
      output <- NULL
      while (length(Lineages)>0){
        i=i+1
        
        scores <- sapply(Lineages,compute_score)
        winner <- which.max(scores)
        
        output_table <- rc_tbl[rc_index %in% Lineages[[winner]]]
        output_table[,Lineage:=i]
        output <- rbind(output,output_table)
        Lineages <- filter_winner(winner,Lineages,row.tree,rc_tbl,row.nodemap)
      }
      return(output)
    }
  } else {
    rc_tbl$Lineage <- 1
    return(rc_tbl)
  }
}
