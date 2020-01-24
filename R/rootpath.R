rootpath <- function(node,tree){
  N=length(tree$tip.label)
  nds <- node
  i=1
  while(!(N+1) %in% nds){
    ix=tree$edge[,2]==nds[i]
    if (!any(ix)){
      stop(paste('Could not find find descendant node',nds[i],'in tree$edge'))
    } else {
      nds <- c(nds,tree$edge[ix,1])
      i=i+1
    }
  }
  return(nds)
}
