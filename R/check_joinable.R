#' check if two rc_seqs are joinable
#' @export
#' @param seq1 index of sequence from \code{\link{rc_seqs}}
#' @param seq2 index of sequence from \code{\link{rc_seqs}}
#' @param Seqs Set of sequences from \code{\link{rc_seqs}}
#' @param rc_table made in \code{\link{makeRCtable}}
#' @param r.nodemap \code{makeNodeMap} of row.tree
#' @param c.nodemap \code{makeNodeMap} of col.tree
check_joinable <- function(seq1,seq2,Seqs,rc_table,r.nodemap,c.nodemap){
  
  ### some of these lists are whole lineages
  ### some are alternative sequences of the same lineage
  
  ### e.g. 
  
  ## 1-2-5-6 and 1-2-3-4 are alternative seqs if they are not lineages
  ## they are lineages if the rc's of {3,5} don't share an 
  ## ancestor-descendant relation in the column/row tree.
  ## The best way to check this will be to find the row/col descendants of 3,5
  ## if their intersect is empty, they are compatible.
  
  ## Alternatively, if we have:
  ## 1-2-3-4 and 1-2-5-4, we know these are alternative seqs 
  ## since 3 and 5 are both on rootpath to 4.
  intrsct <- intersect(Seqs[[seq1]],Seqs[[seq2]])
  if (length(intrsct)==0){
    joinable <- FALSE
  } else {
    disagreements1 <- !Seqs[[seq1]]%in%intrsct
    disagreements2 <- !Seqs[[seq2]]%in%intrsct
    ### if we see FALSE TRUE FALSE - a true surrounded by falses, then these are alternatives
    if (all(c(1,-1) %in% diff(disagreements1)) | 
        all(c(1,-1) %in% diff(disagreements2))){ ### 
      joinable <- FALSE
    } else {
      ix1 <- Seqs[[seq1]][min(which(disagreements1))]
      ix2 <- Seqs[[seq2]][min(which(disagreements2))]
      
      
      nds1 <- unlist(rc_table[rc_index==ix1,c('row.node','col.node')])
      nds2 <- unlist(rc_table[rc_index==ix2,c('row.node','col.node')])
      if (nds1['row.node']==nds2['row.node']){
        joinable <- FALSE
      } else {
        RowDesc1 <- c(unlist(getIndexSets(nds1['row.node'],r.nodemap)),nds1['row.node'])
        RowDesc2 <- c(unlist(getIndexSets(nds2['row.node'],r.nodemap)),nds2['row.node'])
        Common_row_descendants <- intersect(RowDesc1,RowDesc2)
        ColDesc1 <- c(unlist(getIndexSets(nds1['col.node'],c.nodemap)),nds1['col.node'])
        ColDesc2 <- c(unlist(getIndexSets(nds2['col.node'],c.nodemap)),nds2['col.node'])
        Common_col_descendants <- intersect(ColDesc1,ColDesc2)
        n_common_row <- length(Common_row_descendants)
        n_common_col <- length(Common_col_descendants)
        if (n_common_row>0 & n_common_col>0){ ## descendants overlap ==> alternatives
          joinable <- FALSE
        } else { ## descendants don't overlap ==> compatible
          joinable <- TRUE
        }
      }
      
    }
  }
  return(joinable)
}