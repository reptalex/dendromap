#' Simulate multiple treeMaps on row.tree
#' @export
#' @param n max number of independent lineages on \code{row.tree} to map to column tree
#' @param row.tree \code{phylo} class object
#' @param col.tree \code{phylo} class object
#' @param prob.row probability of a node in the row.tree manifesting a map to col.tree. Default is 1
#' @param prob.col probability of a node in the col.tree manifesting a map to row.tree. Default is 1
#' @param sd optional positive numeric, standard deviation input to \code{\link{rnorm}} to create output diagonal matrix
#' @param W optional matrix, \code{\link{treeBasis}} of \code{row.tree}
#' @param V optional matrix, \code{\link{treeBasis}} of \code{col.tree}
#' @param row.node optional start node for row.tree
#' @param col.node optional start node for col.tree
#' @param row.nodeset optional set of nodes to consider for random start of \code{treeMap}
#' @param col.nodeset optional set of nodes to consider for random start of \code{treeMap}
#' @param row.depth.min optional numeric minimum depth for nodes for random start of \code{treeMap}
#' @param row.depth.max optional numeric maximum depth for nodes for random start of \code{treeMap}
#' @param col.depth.min optional numeric minimum depth for nodes for random start of \code{treeMap}
#' @param col.depth.max optional numeric maximum depth for nodes for random start of \code{treeMap}
#' @param fix.col.node logical: whether or not to fix the column node for the start of every \code{treeMap}
#' @param decay.rate exponential rate of decay for past nodes. Probability of draw will be exp(-decay.rate*dT) where dT is the time from column node to row node. Only implemented if \code{use.depths=TRUE}
#' @param use.depth logical: whether or not to use propensities to limit manifestation in row tree to being after nodes in column tree. 
#' @param col.nb \code{\link{nodeBank}} for column tree, allowing user to set up propensities of column tree nodes.
#' @examples 
#' library(ggplot2)
#' library(ggpubr)
#' set.seed(3)
#' m=1e3
#' n=30
#' row.tree <- rtree(m) %>% phytools::force.ultrametric()
#' col.tree <- rtree(n)
#' S <- treeSim(10,row.tree,col.tree,prob.row=0.7,prob.col=0.8,sd = 1e3,
#'              row.depth.min=2,row.depth.max=3)              
#' 
#' plot(S,col.tr.left = 0.47,
#'      col.tr.width = 0.505,
#'      col.tr.bottom = 0.74)
#'  
#' set.seed(1)
#' col.nb <- nodeBank(node=n+1,col.tree,propensity=50,include.node=T)
#' col.nb[,propensity:=50*(max(depth)-depth)]  ##root nodes have higher propensities
#' S <- treeSim(50,row.tree,col.tree,prob.row=0.7,prob.col=0.8,sd = 1e3,
#'              col.nb=col.nb,decay.rate=20,use.depths=T) 
#'              
#' dendroPlot <- plot(S,col.tr.left = 0.47,
#'                    col.tr.width = 0.505,
#'                    col.tr.bottom = 0.74)   
#'      
#' row.depths <- data.table('row.depth'=node.depth.edgelength(row.tree))
#' row.depths[,row.node:=1:.N]
#' setkey(row.depths,row.node)
#' col.depths <- data.table('col.depth'=node.depth.edgelength(col.tree))
#' col.depths[,col.node:=1:.N]
#' setkey(col.depths,col.node)
#' setkey(S$Lineages,row.node)
#' Paths <- row.depths[S$Lineages]
#' setkey(Paths,col.node)
#' Paths <- col.depths[Paths]
#' Paths[,Lineage:=factor(Lineage)]
#' 
#' depthPlot=ggplot(Paths,aes(row.depth,col.depth,color=Lineage))+
#'                  geom_point()+
#'                  geom_abline(intercept=0,slope=1)+
#'                  geom_smooth(method='glm',se=F)+
#'                  scale_x_continuous('Time on Microbial Tree')+
#'                  scale_y_continuous('Time on Reference Tree')
#'           
#' ggarrange(dendroPlot,depthPlot)

treeSim <- function(n,row.tree,col.tree,
                    prob.row=1,prob.col=1,sd=1,W=NULL,V=NULL,
                    row.node=NULL,col.node=NULL,
                    row.nodeset=NULL,col.nodeset=NULL,
                    row.depth.min=NULL,row.depth.max=NULL,
                    col.depth.min=NULL,col.depth.max=NULL,fix.col.node=FALSE,
                    decay.rate=0,use.depths=F,col.nb=NULL,row.nb=NULL){
  
  if (fix.col.node==TRUE & is.null(col.node)){
    col.node <- 1+ape::Ntip(col.tree)
  }
  
  if (is.null(row.depth.min)){
    row.depth.min <- 0
  }
  if (is.null(row.depth.max)){
    row.depth.max <- Inf
  }
  if (is.null(col.depth.min)){
    col.depth.min <- 0
  }
  if (is.null(col.depth.max)){
    col.depth.max <- Inf
  }
  
  depths <- ape::node.depth.edgelength(row.tree)
  if (is.null(row.nodeset)){
    row.nodeset <- intersect(1:ape::Nnode(row.tree)+length(row.tree$tip.label),
                             which(depths<=row.depth.max & depths>=row.depth.min))
  } else {
    row.nodeset <- intersect(row.nodeset,
                             which(depths<=row.depth.max & depths>=row.depth.min))
  }
  
  if (is.null(col.nodeset)){
    col.nodeset <- intersect(1:ape::Nnode(col.tree)+length(col.tree$tip.label),
                             which(depths<=col.depth.max & depths>=col.depth.min))
  } else {
    col.nodeset <- intersect(col.nodeset,
                             which(depths<=col.depth.max & depths>=col.depth.min))
  }
  
  if (is.null(col.nb) & use.depths){
    col.nb=nodeBank(ape::Ntip(col.tree)+1,col.tree,prob.col,include.node=T)
  }
  if (is.null(row.nb) & use.depths){
    row.nb=nodeBank(ape::Ntip(row.tree)+1,row.tree,prob.row,include.node=T)
  }
  
  Path <- vector(mode='list',length=n)
  for (ii in 1:n){
    if (ii==1){
      
      if (is.null(col.node)){
        if (!use.depths){
          col.node <- sample(col.nodeset,1)
        } else {
          lambdas <- col.nb[match(col.nodeset,node),propensity]
          if (any(is.na(lambdas))){
            stop(paste('Erorr in col.nb[match(col.nodeset,node)]: ',sum(is.na(lambdas)),' node(s) from nodeset not found in col.nb',sep=''))
          }
          col.node <- sample(col.nodeset,1,prob=)
        }
      }
      if (is.null(row.node)){
        if (!use.depths){
          row.node <- sample(row.nodeset,1)
        } else {
          row.nb[is.infinite(propensity) & propensity>0,propensity:=max(propensity[!is.infinite(propensity)])]
          if (all(is.infinite(row.nb$propensity))){
            lambdas=rep(1,length(row.nodeset))
          } else {
            lambdas=sapply(row.nb[match(row.nodeset,node),propensity],max,0)
          }
          dT <- row.nb[match(row.nodeset,node),depth]-col.nb[node==col.node,depth]
          lambdas[dT<0]=0
          row.node <- sample(row.nodeset,1,prob=lambdas*exp(-decay.rate*dT))
        }
      }
      
      Path[[ii]] <- treeMap(row.tree,col.tree,row.node,col.node,prob.row,prob.col,decay.rate,col.nb,use.depths)
    } else {
      row.start <- Path[[ii-1]]$row.node[1]
      row.nodeset <- intersect(row.nodeset,disjointNodeset(row.start,row.tree))
      
      if (length(row.nodeset)>0){
        if (!fix.col.node){
          if (!use.depths){
            col.node <- sample(col.nodeset,1)
          } else {
            col.node <- sample(col.nodeset,1,prob=col.nb[match(col.nodeset,node),propensity])
          }
        }
        # row.node <- sample(row.nodeset,1)
        
        
        if (!use.depths){
          row.node <- sample(row.nodeset,1)
        } else {
          row.nb[is.infinite(propensity) & propensity>0,propensity:=max(propensity[!is.infinite(propensity)])]
          if (all(is.infinite(row.nb$propensity))){
            lambdas=rep(1,length(row.nodeset))
          } else {
            lambdas=sapply(row.nb[match(row.nodeset,node),propensity],max,0)
          }
          dT <- row.nb[match(row.nodeset,node),depth]-col.nb[node==col.node,depth]
          lambdas[dT<0]=0
          row.node <- sample(row.nodeset,1,prob=lambdas*exp(-decay.rate*dT))
        }
        Path[[ii]] <-  treeMap(row.tree,col.tree,row.node,col.node,prob.row,prob.col,decay.rate,col.nb,use.depths)
      } else {
        Path <- Path[1:(ii-1)]
        break
      }
    }
  }
  for (i in 1:length(Path)){
    Path[[i]][,Lineage:=i]
  }
  output <- NULL
  output$Lineages <- data.table::rbindlist(Path)
  if (is.null(W)){
    W <- treeBasis(row.tree)
  }
  if (is.null(V)){
    V <- treeBasis(col.tree)
  }
  
  output$W <- W[,output$Lineages$row.node-ape::Ntip(row.tree)]
  dd <- rnorm(nrow(output$Lineages),sd=sd)
  dd <- abs(dd)*output$Lineages$orientation
  output$D <- diag(dd)
  output$V <- V[,output$Lineages$col.node-ape::Ntip(col.tree)]
  output$row.tree <- row.tree
  output$col.tree <- col.tree
  class(output) <- 'dendromap'
  return(output)
}
