#' check if two rc_seqs are joinable
#' @export
#' @param seq1 index of sequence from \code{\link{rc_seqs}}
#' @param seq2 index of sequence from \code{\link{rc_seqs}}
#' @param Seqs Set of sequences from \code{\link{rc_seqs}}
#' @param rc_table made in \code{\link{makeRCtable}}
#' @param Row_Descendants named \code{getIndexSets} of all row.nodes
#' @param Col_Descendants named \code{getIndexSets} of all col.nodes
check_joinable <- function(row.node1,row.node2,col.node1,col.node2,Row_Descendants,Col_Descendants){
  if (row.node1==row.node2 | col.node1==col.node2){
    joinable <- FALSE
  } else {
    RowDesc1 <- c(unlist(Row_Descendants[[as.character(row.node1)]]),row.node1)
    RowDesc2 <- c(unlist(Row_Descendants[[as.character(row.node2)]]),row.node2)
    ColDesc1 <- c(unlist(Col_Descendants[[as.character(col.node1)]]),col.node1)
    ColDesc2 <- c(unlist(Col_Descendants[[as.character(col.node2)]]),col.node2)
    
    Common_row_descendants <- intersect(RowDesc1,RowDesc2)
    Common_col_descendants <- intersect(ColDesc1,ColDesc2)
    n_common_row <- length(Common_row_descendants)
    n_common_col <- length(Common_col_descendants)
    if (n_common_row>0 & n_common_col>0){ ## descendants overlap ==> altenatives
      joinable <- FALSE
    } else { ## descendants don't overlap ==> compatible
      joinable <- TRUE
    }
  }
  return(joinable)
}
