#' grow_lineage from row/col edge pair. internal function for rdendromap
#' @export
#' @param row.edge
#' @param col.edge
#' @param row.tree
#' @param col.tree
#' @param rowEdgeMap
#' @param colEdgeMap
#' @param lineage_size
#' @param rowDescendants
#' @param colDescendants
grow_lineage <- function(row.edge,col.edge,row.tree,col.tree,
                         rowEdgeMap,colEdgeMap,lineage_size,
                         rowDescendants,colDescendants){
  output <- data.table('row.edge'=row.edge,'col.edge'=col.edge)
  row.edges <- rowDescendants[[row.edge]]
  col.edges <- colDescendants[[col.edge]]
  if (lineage_size==1 | length(row.edges)==0){
    done=TRUE
  } else {
    done=FALSE
  }
  i=1
  

  rc_tbl <- expand.grid('row.edge'=row.edges,'col.edge'=col.edges) %>% as.data.table
  rc_tbl[,rc_index:=1:.N]
  while(i<lineage_size & nrow(rc_tbl)>0){
    
    ix <- sample(nrow(rc_tbl),1)
    
    re <- rc_tbl$row.edge[ix]
    ce <- rc_tbl$col.edge[ix]
    output <- rbind(output,data.table('row.edge'=re,'col.edge'=ce))
    
    rea <- edge_ancestors(re,row.tree)
    cea <- edge_ancestors(ce,col.tree)
    
    #### remove incompatibles #### 
    rc_ix <- rc_tbl$rc_index[ix]
    row.edg <- rc_tbl[rc_index==rc_ix,row.edge]
    col.edg <- rc_tbl[rc_index==rc_ix,col.edge]
    row_descs <- rowDescendants[[row.edg]]
    col_descs <- colDescendants[[col.edg]]
    row_ancs <- c(edge_ancestors(row.edg,row.tree),row.edg)
    incompatibles <- rc_tbl[row.edge %in% row_ancs,rc_index]
    if (is.null(col_descs)){ ## No possible descendants
      incompatibles <- c(incompatibles,rc_tbl[row.edge %in% row_descs,rc_index])
    } else {
      incompatibles <- c(incompatibles,rc_tbl[row.edge %in% row_descs & !col.edge %in% col_descs,rc_index])
    }
    rc_tbl <- rc_tbl[!rc_index %in% incompatibles]
    #############################
    i=i+1
  }
  return(output)
}

