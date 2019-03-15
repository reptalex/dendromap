manifestRowChildren <- function(node,row.tree,row.nb,col.nb,col.nodes,use.depths=F){
  if (use.depths){
    colnds <- c(col.nodes,unlist(phangorn::Descendants(col.tree,col.nodes,'all')))
    colnds <- colnds[colnds>ape::Ntip(col.tree)]
    col.depths <- col.nb[node%in%colnds,c('node','depth','propensity')]
    NB <- row.nb[depth>=min(col.depths)]
    
    if (nrow(NB)==0){
      MAP <- NULL
    } else {
      if (node %in% NB$node){
        m <- node
      } else {
        snb <- row.nb
        snb[,manifest:=FALSE]
        snb[node %in% NB$node,manifest:=TRUE]
        m <- manifestChildren(node,row.tree,snb)
      }
      MAP <- NULL
      while (length(m)>0){
        
        map <- depthDraw(NB[match(m,node),depth],col.depths$depth,col.depths$propensity)
        map[,row.node:=m[row]]
        m <- map[col==0,row.node]
        map <- map[col>0]
        map[,col.node:=col.depths$node[col]]
        map[,orientation:=sign(rnorm(.N))]
        map[,terminated:=row.node %in% row.nb[terminal==T,node] | 
                         col.node %in% col.nb[terminal==T,node]]
        MAP <- rbind(MAP,map[,c('row.node','col.node','orientation','terminated')])
        if (length(m)>0){
          ch <- sapply(m,labelledChildren,row.tree)
          m <- c(ch)
          m <- m[m>ape::Ntip(row.tree)]
        }
      }
    }
  } else {
    
    MAP <- data.table('row.node'=manifestChildren(nd = node,tree=row.tree,nodebank=row.nb))
    if (length(col.nodes)>1){
      MAP[,col.node:=sample(col.nodes,nrow(MAP),T,prob=col.nb[match(col.nodes,node),propensity])]
    } else {
      MAP[,col.node:=col.nodes]
    }
    MAP[,orientation:=sign(rnorm(.N))]
    MAP[,terminated:=row.node %in% row.nb[terminal==T,node] | 
          col.node %in% col.nb[terminal==T,node]]
  }
  return(MAP)
}
