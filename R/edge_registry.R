#' Make registry of edges and number of descendants for quick traversals
#' @export
#' @param tree phylo class object
#' @examples
edge_registry <- function(tree){
  ## while nodemap had pos/neg, this will have edge,n = number descendants
  ntip <- ape::Ntip(tree)
  dt <- data.table('edge'=1:nrow(tree$edge),'n'=0)
  setkey(dt,edge)
  
  root_edges <- which(tree$edge[,1]==(ape::Ntip(tree)+1))
  
  tips <- data.table('edge'=which(tree$edge[,2]<=length(tree$tip.label)))
  tips[,anc:=match(tree$edge[edge,1],tree$edge[,2])]
  tips <- tips[!is.na(anc)] ## this is a scenario where a root branch is a tip.
  ## we'll loop through tips, joining them, defining their ancestor as tip.
  while(nrow(tips)>0){
    common_ancs <- as.numeric(names(which(table(tips$anc)==2)))
    jt <- tips[anc %in% common_ancs,list(edge=unique(anc),
                                         d1=edge[1],
                                         d2=edge[2]),by=anc]
    
    dt[match(jt$edge,edge),n:=n+dt[match(jt$d1,edge),n]+dt[match(jt$d2,edge),n]+2]
    tips <- tips[! anc %in% common_ancs]
    new_tips <- data.table('edge'=setdiff(common_ancs,root_edges))
    new_tips[,anc:=match(tree$edge[edge,1],tree$edge[,2])]
    tips <- rbind(tips,new_tips)
  }
  setkey(dt,edge)
  return(dt)
}
