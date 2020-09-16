library(phylofactor)
library(dendromap)
library(parallel)


edge_rc_table <- function(X,row.tree,col.tree,maxPval=0.0001,col_shuffle_only=FALSE){
  m <- nrow(X)
  n <- ncol(X)
  W <- getPhyloGroups(row.tree) %>% sapply(ilrvec,n=nrow(X))
  V <- getPhyloGroups(col.tree) %>% sapply(ilrvec,n=ncol(X))
  
  U <- t(W) %*% X %*% V
  if (col_shuffle_only){
    Unull <- t(W) %*% X[,sample(n)] %*% V
  } else {
    Unull <- t(W) %*% X[sample(m),sample(n)] %*% V
  }
  
  cdf <- ecdf(c(Unull))
  
  
  rc_table=which(abs(U)>0) %>% arrayInd(.dim = c(nrow(U),ncol(U))) %>% as.data.table
  names(rc_table) <- c('row.edge','col.edge')
  rc_table$stat <- U[abs(U)>0]
  setkey(rc_table,row.edge,col.edge,stat)
  
  rc_table[,rank:=rank(-abs(stat))]
  rc_table[,P:=1-cdf(stat)]
  rc_table <- rc_table[stat>0]
  setkey(rc_table,stat)
  
  y <- log(c(Unull[Unull>0]))
  if (any(rc_table$P==0)){
    nn=sum(rc_table$P==0)
    base::cat(paste('\n',nn,' P-values were 0. Will estimate tail probabilities assuming log(stat^2)~rnorm(mu,sd)',sep=''))
    y <- y[y>-20]
    mu <- mean(y)
    sig <- sd(y)
    min.P <- rc_table[P>0,min(P)]
    stat.min <- max(log(rc_table[P>0][P==min(P),stat]))
    estimate_tail <- function(y1,mu,sig,ymin,pmin) (1-pnorm(y1,mu,sig))/(1-pnorm(ymin,mu,sig))*pmin
    rc_table[P==0,P:=estimate_tail(log(stat),mu,sig,stat.min,min.P)]
  }
  
  rc_table <- rc_table[P<=maxPval]
  rc_table[,rc_index:=1:.N]
  
  return(rc_table)
}

makeEdgeMap <- function(tree){
  ## while nodemap had pos/neg, this will have edge,n = number descendants
  ntip <- ape::Ntip(tree)
  dt <- data.table('edge'=1:nrow(tree$edge),'n'=0)
  setkey(dt,edge)
  
  root_edges <- which(tree$edge[,1]==(ape::Ntip(tree)+1))
  
  tips <- data.table('edge'=which(tree$edge[,2]<=length(tree$tip.label)))
  tips[,anc:=match(tree$edge[edge,1],tree$edge[,2])]
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

edgeRCMap <- function(rc_table,Row_Descendants,Col_Descendants){
  
  row.terminals <- which(sapply(Row_Descendants,FUN=function(x) length(x)==0))
  col.terminals <- which(sapply(Col_Descendants,FUN=function(x) length(x)==0))
  
  ix <- which((!rc_table$row.edge %in% row.terminals) &
                (!rc_table$col.edge %in% col.terminals))
  maps <- lapply(ix,makeDescendantTable,
                 rc_table,
                 row.terminals,
                 col.terminals,
                 Row_Descendants,
                 Col_Descendants,
                 method='edge')
  RCmap <- rbindlist(maps)
  RCmap[,terminal:=!(descendant %in% ancestor)]
  setkey(RCmap,descendant,ancestor)
  return(RCmap)
}

edge_ancestors <- function(edge,tree){
  rt=ape::Ntip(tree)+1
  anc=tree$edge[edge,1]
  if (anc==rt){
    return(NULL)
  } else {
    anc_edges <- NULL
    while(!anc==rt){
      new_anc_edge <- which(tree$edge[,2]==anc)
      anc_edges <- c(anc_edges,new_anc_edge)
      anc=tree$edge[new_anc_edge,1]
    }
    return(anc_edges)
  }
}

######## option 1: greedyColTree
incompatible_descendants <- function(rc_ix,rc_tbl,Row_Descendants,Col_Descendants,row.tree.=row.tree){
  row.edg <- rc_tbl[rc_index==rc_ix,row.edge]
  col.edg <- rc_tbl[rc_index==rc_ix,col.edge]
  ### incompatible lineages will have either
  ## the same col.edges in descendant row.edges
  ## or ancestral edges
  
  row_descs <- Row_Descendants[[row.edg]]
  col_descs <- Col_Descendants[[col.edg]]
  
  row_ancs <- edge_ancestors(row.edg,row.tree)
  incompatibles <- rc_tbl[(col.edge==col.edg & row.edge %in% c(row_descs,row_ancs))|
                              row.edge %in% row_descs & !col.edge %in% col_descs,rc_index]
  incompatibles <- c(incompatibles,rc_table[row.edge==row.edg,rc_index])
  return(incompatibles)
}

clean_sisters <- function(lineage,row.tree,rc_table){
  ### any sister row.edges with same col.edge will reclassify their mother as their col.edge
  repeated_col_edges <- as.numeric(names(which(table(lineage$col.edge)>0)))
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

greedyColTree <- function(rc_tbl,RCmap.=RCmap,rc_table.=rc_table,row.tree.=row.tree,method='basal'){
  lineage=NULL
  done=F
  n=0
  setkey(rc_tbl,row.edge,col.edge)
  while(!done){
    n=n+1
    lineage <- rbind(lineage,rc_tbl[1,])
    if (method=='basal'){
      ix <- rc_tbl$rc_index[1]
    } else {
      ix <- rc_tbl[P==min(P),rc_index[1]]
    }
    if (n==1){
      setkey(rc_tbl,col.edge,P)
    }
    incompatibles <- c(ix,incompatible_descendants(ix,rc_tbl,Row_Descendants,Col_Descendants))
    rc_tbl <- rc_tbl[!rc_index %in% incompatibles]
    if (nrow(rc_tbl)==0){
      done=TRUE
    }
  }
  lineage <- clean_sisters(lineage,row.tree,rc_table)
  lineage[,lineage_id:=rc_index[1]]
  return(lineage)
}

getLineage <- function(basal_ix,RCmap.=RCmap,rc_table.=rc_table,row.tree.=row.tree,
                       Row_Descendants.=Row_Descendants,Col_Descendants.=Col_Descendants){
  basal_row_edge <- rc_table[rc_index==basal_ix,row.edge]
  basal_col_edge <- rc_table[rc_index==basal_ix,col.edge]
  row_descs <- Row_Descendants[[basal_row_edge]]
  col_descs <- Col_Descendants[[basal_col_edge]]
  rc_tbl <- rc_table[row.edge %in% c(basal_row_edge,row_descs) & 
                       col.edge %in% c(basal_col_edge,col_descs)]
  lineage <- greedyColTree(rc_tbl,RCmap,rc_table,row.tree)
  lineage[,lineage_id:=basal_ix]
  setkey(lineage,col.edge,P)
  return(lineage)
}


edgeTips <- function(tree){
  grps <- phylofactor::getPhyloGroups(tree)
  mat <- lapply(grps,FUN=function(gg)
    data.table('min'=min(gg[[1]]),'max'=max(gg[[1]])))
  mat <- rbindlist(mat)
  mat[,edge:=1:.N]
  setkey(mat,edge)
  mat[,tip:=(min==max)]
  return(mat)
}

lineageBoxes <- function(lineage,X,colEdgeTips,rowEdgeTips){
  #### if we start at the finest lineage w/ highest row.edge, we can construct bins iteratively
  #### 
  setkey(lineage,row.edge,col.edge)
  rowset <- rowEdgeTips[edge==lineage$row.edge[1],seq(min,max)]
  colset <- colEdgeTips[edge==lineage$col.edge[1],seq(min,max)]
  
  Xdt <- data.table('i'=rep(rowset,times=ncol(X)),
                    'j'=rep(1:ncol(X),each=length(rowset)),
                    'x'=c(X[rowset,]),
                    'lineage_id'=unique(lineage$lineage_id),
                    'box'=0)
  
  setkey(lineage,row.edge)
  setkey(Xdt,i,j)
  
  #### rowset will be continuously whittled down
  #### colset will need to be iterated.
  mbox=0
  for (nn in 1:nrow(lineage)){
    mbox=mbox+1
    if (nn==1){
      Xdt[j %in% colset,box:=mbox]
    } else {
      spp <- rowEdgeTips[edge==lineage$row.edge[nn],seq(min,max)]
      cols <- colEdgeTips[edge==lineage$col.edge[nn],seq(min,max)]
      bx <- Xdt[i %in% spp & box!=0,unique(box)] # The first box is the full set of species in the continents in which they're not found
      
      Xdt[i %in% spp & box!=0,box:=mbox]
      mbox=mbox+1
      Xdt[i %in% spp & j %in% cols,box:=mbox]
    }
  }
  return(Xdt)
}

######## THIS TAKES A WHILE - PARALLELIZE 
lineage_stats <- function(Lineages,X,colEdgeTips,rowEdgeTips){
  lineage_ids <- unique(Lineages$lineage_id)
  fits <- data.table('lineage_id'=lineage_ids,
                     'F_stat'=0)
  for (id in lineage_ids){
    Xdt <- lineageBoxes(Lineages[lineage_id==id],X,colEdgeTips,rowEdgeTips)
    setkey(Xdt,box)
    Xdt[,box:=factor(box)]
    a <- aov(x~box,data=Xdt)
    ss <-summary(a)[[1]]['box',] %>% unlist
    fits[lineage_id==id,F_stat:=ss['F value']]
    fits[lineage_id==id,msq:=ss['Mean Sq']]
    fits[lineage_id==id,ssq:=ss['Sum Sq']]
  }
  return(fits)
}

global_Fstat <- function(Lineages,X,colEdgeTips,rowEdgeTips){
  lineages <- split(Lineages,f = factor(Lineages$lineage_id)) %>% 
    lapply(lineageBoxes,X,colEdgeTips,rowEdgeTips) %>% rbindlist
  
  Xdt <- data.table('i'=rep(1:nrow(X),times=ncol(X)),
                    'j'=rep(1:ncol(X),each=nrow(X)),
                    'x'=c(X))
  setkey(Xdt,i,j)
  setkey(lineages,i,j)
  Xdt <- lineages[,c('i','j','lineage_id','box')][Xdt]
  Xdt[is.na(lineage_id),lineage_id:=-1]
  Xdt[is.na(box),box:=-1]
  Xdt[,group:=paste(lineage_id,box,sep='_')]
  Xdt[,yhat:=mean(x,na.rm=T),by=group]
  a=aov(x~box,data=Xdt) %>% summary
  return(a[[1]]['box','F value'])
}

filter_stats <- function(x,Lineages,rc_table,RCmap,Row_Descendants){
  x <- x[order(F_stat,decreasing = T)]
  
  i=1
  while(i<=nrow(x)){
    id <- x$lineage_id[i]
    row.edg <- rc_table[rc_index==id,row.edge]
    descendants <- Row_Descendants[[row.edg]]
    ancestors <- edge_ancestors(row.edg,row.tree)
    incompatible_lineages <- setdiff(Lineages[row.edge %in% c(descendants,ancestors),lineage_id],id)
    x <- x[!lineage_id %in% incompatible_lineages]
    i=i+1
  }
  return(x)
}

edge2node <- function(edges,tree) tree$edge[edges,2]

predict.dendromap <- function(object,...){
  
  lineages <- split(object$Lineages,f = factor(object$Lineages$lineage_id)) %>% 
    lapply(lineageBoxes,object$Data,object$colEdgeTips,object$rowEdgeTips) %>% rbindlist
  
  Xdt <- data.table('i'=rep(1:nrow(X),times=ncol(X)),
                    'j'=rep(1:ncol(X),each=nrow(X)),
                    'x'=c(X))
  setkey(Xdt,i,j)
  setkey(lineages,i,j)
  Xdt <- lineages[,c('i','j','lineage_id','box')][Xdt]
  Xdt[is.na(lineage_id),lineage_id:=-1]
  Xdt[is.na(box),box:=-1]
  Xdt[,group:=paste(lineage_id,box,sep='_')]
  Xdt[,yhat:=mean(x,na.rm=T),by=group]
  Xdt <- matrix(Xdt$yhat,nrow=nrow(X),ncol=ncol(X),byrow=T)
  rownames(Xdt) <- rownames(object$Data)
  colnames(Xdt) <- colnames(object$Data)
  return(Xdt)
  
}

dendromap.e <- function(X,row.tree,col.tree,maxPval=0.0001){
  # edge rc_table -----------------------------------------------------------
  rc_table <- edge_rc_table(X,row.tree,col.tree,maxPval)
  
  # edgeMap for quick descendant calculation -------------------------------------------------------------------
  rowEdgeMap <- makeEdgeMap(row.tree)
  colEdgeMap <- makeEdgeMap(col.tree)
  
  
  # RC_map ------------------------------------------------------------------
  row.edges <- 1:Nedge(row.tree)
  col.edges <- 1:Nedge(col.tree)
  Row_Descendants <- rowEdgeMap[,list(list(setdiff(seq(edge,edge+n),edge))),by=edge]$V1
  Col_Descendants <- colEdgeMap[,list(list(setdiff(seq(edge,edge+n),edge))),by=edge]$V1
  RCmap <- edgeRCMap(rc_table,Row_Descendants,Col_Descendants)
  
  # choosing edges to input to greedyColTree -----------------------------------------
  
  basal_indexes <- RCmap[,setdiff(unique(ancestor),unique(descendant))]
  desc_count <- RCmap[,list(n=.N),by=ancestor]
  basal_ixs_with_descendants <- intersect(basal_indexes,desc_count[n>1,ancestor])
  
  cl <- parallel::makeCluster(7)
  
  parallel::clusterExport(cl,varlist=c('incompatible_descendants','clean_sisters','greedyColTree','getLineage','edge_ancestors',
                                       'RCmap','rc_table','Row_Descendants','Col_Descendants','row.tree','col.tree'),
                          envir = environment())
  parallel::clusterEvalQ(cl,library(data.table))
  parallel::clusterEvalQ(cl,library(magrittr))
  Lineages <- parallel::parLapply(cl,sample(basal_ixs_with_descendants),getLineage) %>% rbindlist
  parallel::stopCluster(cl)
  rm('cl')
  
  Lineages[,n:=.N,by=lineage_id]
  Lineages <- Lineages[n>1]
  
  
  # Filtering lineages ------------------------------------------------------
  
  colEdgeTips <- edgeTips(col.tree)
  rowEdgeTips <- edgeTips(row.tree)
  
  ### get fit stats
  stats <- lineage_stats(Lineages,X,colEdgeTips,rowEdgeTips)
  stats <- stats[order(F_stat,decreasing = T)]
  stats <-  filter_stats(stats,Lineages,rc_table,RCmap,Row_Descendants)
  
  Lineages <- Lineages[lineage_id %in% stats$lineage_id]
  
  setkey(Lineages,lineage_id,row.edge)
  setkey(stats,lineage_id)
  
  Lineages <- stats[,c('lineage_id','F_stat')][Lineages]
  setkey(Lineages,row.edge)
  Lineages[,row.node:=edge2node(row.edge,row.tree)]
  Lineages[,col.node:=edge2node(col.edge,col.tree)]
  Lineages <- Lineages[order(F_stat,decreasing = T)]
  Lineages[,Lineage:=match(lineage_id,unique(lineage_id))]
  object <- list('Lineages'=Lineages,
                 'Data'=X,
                 'row.tree'=row.tree,
                 'col.tree'=col.tree,
                 'colEdgeTips'=colEdgeTips,
                 'rowEdgeTips'=rowEdgeTips,
                 'RowDescendants'=Row_Descendants,
                 'ColDescendants'=Col_Descendants)
  class(object) <- 'dendromap'
  return(object)
}

plot.dendromap <- function(x,y=NULL,
                           color.fcn.clade=viridis::viridis,
                           color.fcn.node=viridis::viridis,
                           heatmap.offset=0,
                           col.tr.left=0.5,
                           col.tr.width=0.45,
                           col.tr.bottom=0.75){
  
  # if (is.null(x$Lineages$orientation)){
  #   x$Lineages[,orientation:=sign(stat)]
  # }
  if (any(is.na(x$row.tree$edge.length))){
    x$row.tree$edge.length[is.na(x$row.tree$edge.length)] <- 0
  }
  if (any(is.na(x$col.tree$edge.length))){
    x$col.tree$edge.length[is.na(x$col.tree$edge.length)] <- 0
  }
  
  
  
  vcols <- color.fcn.node(length(unique(x$Lineages$col.node)))
  nodecols <- data.table('node'=sort(unique(x$Lineages$col.node),decreasing = F),
                         'color'=vcols)
  gtr <- ggtree::ggtree(x$col.tree,branch.length = 'none')+
    ggtree::geom_point2(ggplot2::aes(subset=node %in% nodecols$node),
                        color=nodecols$color,cex=3)+
    ggplot2::coord_flip()+ggplot2::scale_x_reverse()+
    ggplot2::scale_y_reverse()
  
  ###### processing data matrix
  column.order <- gtr$data$label[order(gtr$data$y[1:ape::Ntip(x$col.tree)],decreasing = T)]
  if (is.null(y)){
    if (is.null(x$Data)){
      if (is.null(x$W) | is.null(x$V) | is.null(x$D)){
        y <- base::matrix(0,nrow=ape::Ntip(x$row.tree),ncol=ape::Ntip(col.tree))
      } else {
        y <- x$W %*% x$D %*% t(x$V)
        y[y==0] <- NA
      }
    } else {
      y <- x$Data
    }
  }
  
  gg <- ggtree::ggtree(x$row.tree,layout='rectangular',branch.length = 'none')
  
  ## add clade hilights for basal nodes
  cols <- color.fcn.clade(length(unique(x$Lineages$lineage_id)))
  ii=0
  start.nodes <- x$Lineages[,list(nd=min(row.node)),by=lineage_id]$nd
  for (nd in start.nodes){
    ii=ii+1
    gg <- gg+ggtree::geom_hilight(nd,fill=cols[ii])
  }
  
  Path <- x$Lineages
  setkey(Path,row.node)
  gg.cols <- nodecols[match(x$Lineages$col.node,node),]$color
  gg <- gg+ggtree::geom_point2(ggplot2::aes(subset=node %in% Path$row.node),
                               color=gg.cols,cex=3)
  
  gg <- ggtree::gheatmap(gg,as.matrix(y),color=NA,colnames = FALSE,offset=heatmap.offset,high='steelblue',low='grey')+
    ggtree::theme(legend.position = 'none')+cowplot::theme_nothing()
  
  output <- cowplot::ggdraw() + 
    cowplot::draw_plot(gtr, x = col.tr.left, y = col.tr.bottom,
                       width = col.tr.width, height = .25)+
    cowplot::draw_plot(gg, x = 0, y = 0.05, width = 1, height = .7)
  
  return(output)
}


lineage_plot <- function(dm,lineage=NULL,id=NULL,...){
  ### trim row.tree to only lineage; rename/follow edges for right labelling
  ### plot subset of data
  
  if (is.null(id)){
    L <- dm$Lineages[Lineage==lineage]
  } else {
    L <- dm$Lineages[lineage_id==id]
  }
  
  
  basal_edge <- L[,min(row.edge)]
  desc_edges <- dm$RowDescendants[[basal_edge]]
  tips <- phangorn::Descendants(dm$row.tree,dm$row.tree$edge[basal_edge,2],'tips')[[1]]
  rt <- ape::drop.tip(dm$row.tree,setdiff(dm$row.tree$tip.label,dm$row.tree$tip.label[tips]))
  rownames(rt$edge) <- desc_edges
  
  L[row.edge!=basal_edge,row.node:=rt$edge[as.character(row.edge),2]]
  L[row.edge==basal_edge,row.node:=(length(rt$tip.label)+1)]
  
  temp_dm <- list('Data'=dm$Data[rt$tip.label[tips],],
                  'Lineages'=L,
                  'col.tree'=dm$col.tree,
                  'row.tree'=rt)
  class(temp_dm) <- 'dendromap'
  pl <- plot.dendromap(temp_dm,...)
  return(pl)
}
