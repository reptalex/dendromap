#' Walk ant path from descendant - get sequence of ancestral rc's from distal rc
#' @export
#' @param i index - wihch descendant of \code{RCmap} to for which to find ancestral rc sequence
#' @param RCmap made in \code{\link{makeRCmap}}

rc_seqs <- function(i,RCmap){
  nd <- RCmap$descendant[i]
  
  ### note: a given node may have multiple paths towards the root
  ### e.g. it's possible to have nodes (1,2,3,4,nd)
  ### such that either (1,3,nd) or (2,4,nd) are sequences due to
  ### incompatible orientations of 1-2, 3-4, etc.
  ### Hence, for a given terminal node, we may need to make a list
  ### of possible sequences.
  ### junctures will be identifiable based on 4-nd and 3-nd relations
  ### but no 3-4 or 4-3 relation. At such a juncture, we'll replicate
  ### the rc.seq currently rooted at 'nd' and then walk up one juncture,
  ### and when we're done we'll return to walk up the other branch.
  ### We'll have to keep track of the index on node.seq and its un-walked branch on RCmap
  
  ancs <- sort(RCmap[descendant==nd,ancestor],decreasing = T)
  if (length(ancs)>1){
    ### Some of these ancs will be sequences
    ### others will be incompatible
    ### e.g. if we have 1-nd and 2-nd, 
    ### and 1-2 is in RCmap, then we know P(1,2,nd) dominates as a longer alignment
    ### thus we don't need to consider 1-nd branching.
    ### recall that lower indexes correspond to lower column nodes
    ### hence the highest index will not have descendants.
    jj=0
    seq <- NULL
    while(length(ancs)>0){
      searching=T
      jj=jj+1
      seq[[jj]] <- nd
      temp_ancs <- ancs
      while(searching==T){
        seq[jj] <- list(c(temp_ancs[1],seq[[jj]]))
        temp_ancs <- sort(RCmap[descendant==temp_ancs[1],ancestor],decreasing = T)
        if (length(temp_ancs)==0){
          searching=F
        }
      }
      ancs <- setdiff(ancs,seq[[jj]])
    }
  } else {
    seq <- list(c(ancs,nd))
  }
  return(seq)
}