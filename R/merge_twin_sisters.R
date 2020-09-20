#' Identify sister row edges with same column edge association & merge at mother
#' @export
#' @param lineage lineage table. See \code{\link{get_lineage}}
#' @param row.tree \code{phylo} class object
#' @param col.tree \code{phylo} class object
#' @examples
merge_twin_sisters <- function(lineage,row.tree,rc_table){
  ### any sister row.edges with same col.edge will reclassify their mother as their col.edge
  repeated_col_edges <- as.numeric(names(which(table(lineage$col.edge)>1)))
  if (length(repeated_col_edges)==0){
    return(lineage)
  } else {
    row_edges <- lineage[col.edge %in% repeated_col_edges,row.edge]
    ancs <- row.tree$edge[row_edges,1]
    while(any(table(ancs)>1)){
      common_ancs <- as.numeric(names(which(table(ancs)==2)))
      for (anc in common_ancs){
        anc_edge <- which(row.tree$edge[,2]==anc)
        desc_edges <- row_edges[ancs==anc]
        col_edge <- lineage[row.edge==desc_edges[1],col.edge]
        lineage <- lineage[! row.edge %in% c(anc_edge,desc_edges)]
        if (nrow(rc_table[row.edge==anc_edge & col.edge==col_edge])==1){
          lineage <- rbind(lineage,rc_table[row.edge==anc_edge & col.edge==col_edge])
        } else {
          lineage <- rbind(lineage,
                           data.table('row.edge'=anc_edge,'col.edge'=col_edge,
                                      'stat'=NA,'rank'=NA,'P'=NA,'rc_index'=NA))
        }
      }
      row_edges <- lineage[col.edge %in% repeated_col_edges,row.edge]
      ancs <- row.tree$edge[row_edges,1]
    }
    return(lineage)
  }
}