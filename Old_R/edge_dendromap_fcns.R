library(phylofactor)
library(dendromap)
library(parallel)
anc_desc_table <- function(i=1,rc_table,
                                terminalRowNodes,
                                terminalColNodes,
                                rowDescendants,
                                colDescendants,
                                method='node'){
  if (method=='node'){
    row.node <- rc_table$row.node[i]
    col.node <- rc_table$col.node[i]
    if (row.node %in% terminalRowNodes | 
        col.node %in% terminalColNodes){
      return(NULL)
    } else {
      sgn <- sign(rc_table$stat[i])
      if (sgn==1){ ## pos-pos and neg-neg
        rowdescPos <- rowDescendants[[toString(row.node)]][['pos']]
        coldescPos <- colDescendants[[toString(col.node)]][['pos']]
        rowdescNeg <- rowDescendants[[toString(row.node)]][['neg']]
        coldescNeg <- colDescendants[[toString(col.node)]][['neg']]
      } else {
        rowdescPos <- rowDescendants[[toString(row.node)]][['pos']]
        coldescPos <- colDescendants[[toString(col.node)]][['neg']]
        rowdescNeg <- rowDescendants[[toString(row.node)]][['neg']]
        coldescNeg <- colDescendants[[toString(col.node)]][['pos']]
      }
      ixs <- rc_table[(row.node %in% rowdescPos & col.node %in% coldescPos) |
                        (row.node %in% rowdescNeg & col.node %in% coldescNeg)]$rc_index
      if (length(ix)>0){
        return(data.table('ancestor'=rc_table$rc_index[i],
                          'descendant'=ixs))
      } else {
        return(NULL)
      }
    }
  } else {
    row.edge <- rc_table$row.edge[i]
    col.edge <- rc_table$col.edge[i]
    if (row.edge %in% terminalRowNodes|
        col.edge %in% terminalColNodes){
      return(NULL)
    } else {
      
      rowDesc <- rowDescendants[[row.edge]]
      colDesc <- colDescendants[[col.edge]]
      
      ixs <- rc_table[(row.edge %in% rowDesc) & (col.edge %in% colDesc)]$rc_index
      if (length(ixs)>0){
        return(data.table('ancestor'=rc_table[i]$rc_index,
                          'descendant'=ixs))
      } else {
        return(NULL)
      }
    }
  }
}

make_rc_table <- function(X,row.tree,col.tree,maxPval=0.0001,col_shuffle_only=FALSE){
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

edge_registry <- function(tree){
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

find_rc_relations <- function(rc_table,rowDescendants,colDescendants){
  
  row.terminals <- which(sapply(rowDescendants,FUN=function(x) length(x)==0))
  col.terminals <- which(sapply(colDescendants,FUN=function(x) length(x)==0))
  
  ix <- which((!rc_table$row.edge %in% row.terminals) &
                (!rc_table$col.edge %in% col.terminals))
  maps <- lapply(ix,anc_desc_table,
                 rc_table,
                 row.terminals,
                 col.terminals,
                 rowDescendants,
                 colDescendants,
                 method='edge')
  rc_relations <- rbindlist(maps)
  rc_relations[,terminal:=!(descendant %in% ancestor)]
  setkey(rc_relations,descendant,ancestor)
  return(rc_relations)
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

incompatible_descendants <- function(rc_ix,rc_tbl,rowDescendants,colDescendants,row.tree.=row.tree,col.tree.=col.tree){
  row.edg <- rc_tbl[rc_index==rc_ix,row.edge]
  col.edg <- rc_tbl[rc_index==rc_ix,col.edge]
  ### incompatible lineages will have either
  ## the same (or ancestral) col.edges in descendant row.edges
  ## or ancestral row.edges
  
  ## in other words, row_ancs are excluded and only descendant col.edges are allowed in descednant row edges
  
  row_descs <- rowDescendants[[row.edg]]
  col_descs <- colDescendants[[col.edg]]
  row_ancs <- c(edge_ancestors(row.edg,row.tree),row.edg)
  incompatibles <- rc_tbl[row.edge %in% row_ancs,rc_index]
  if (is.null(col_descs)){ ## No possible descendants
    incompatibles <- c(incompatibles,rc_tbl[row.edge %in% row_descs,rc_index])
  } else {
    incompatibles <- c(incompatibles,rc_tbl[row.edge %in% row_descs & !col.edge %in% col_descs,rc_index])
  }
  return(incompatibles)
}



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



assemble_lineage <- function(rc_tbl,rc_relations.=rc_relations,rc_table.=rc_table,row.tree.=row.tree,rowDescendants.=rowDescendants,filter=TRUE){
  ### RULES:
  ## 1) remove col.edges with more significant values among descendants (implies signal lost with inclusion of sister taxa)
  ##     see function basal_dominance_filter
  ## 2) pick most basal row.edge, choosing col.edge with lowest P-value  -- from filtering, this is guaranteed to be greater than that of any descendant for same col.edge
  ## 3) remove row.edge and all descendant rc_indexes with same col.edge from rc_tbl
  ## 4) repeat 2-3 until rc_tbl is empty
  
  lineage=NULL
  done=F
  n=0
  setkey(rc_tbl,row.edge,P)
  while(!done){
    n=n+1
    lineage <- rbind(lineage,rc_tbl[1,]) ## next most basal row edge, and its lowest P-value col.edge
    ix <- rc_tbl$rc_index[1]
    incompatibles <- c(ix,incompatible_descendants(ix,rc_tbl,rowDescendants,colDescendants))
    rc_tbl <- rc_tbl[!rc_index %in% incompatibles]
    if (nrow(rc_tbl)==0){
      done=TRUE
    }
  }
  lineage <- merge_twin_sisters(lineage,row.tree,rc_table)
  return(lineage)
}

basal_dominance_filter <- function(rc_tbl,rowDescendants){
  ### Filtering ### 
  ## We'd like to do this as few times as possible.
  ## One option is to start with basal indexes of rc_table for the min P-value cutoff
  ## and then use that filtered rc_table's basal-index shards for assemble_lineage input
  ## When adding more rc-indexes with higher P-values, 
  ## we can simply check whether/not they, too, need to be excluded.
  
  ## RULES OF FILTERING:
  ## remove col.edges with more significant values among descendants (implies signal lost with inclusion of sister taxa)
  ### Implications:
  ## addition of P-values with higher values will not upset a lower P-value's filtering. 
  setkey(rc_tbl,row.edge,P)
  row.edges <- rc_tbl[,unique(row.edge)]
  for (ee in row.edges){
    dum <- rc_tbl[row.edge==ee,c('col.edge','P','rc_index')]
    if (nrow(dum)>0 & length(rowDescendants[[ee]])>0){
      descs_tbl <- rc_tbl[row.edge %in% rowDescendants[[ee]] & col.edge %in% dum$col.edge]
      descs <- descs_tbl[,list(Pmin=min(P)),by=col.edge]
      setkey(dum,col.edge)
      setkey(descs,col.edge) ## this may not contain all col.edges in dum
      dd <- descs[dum]
      removers <- dd[P>Pmin,rc_index]
      rc_tbl <- rc_tbl[!rc_index %in% removers]
    }
  }
  return(rc_tbl)
}
simplify_rc_table <- function(rc_table,rc_relations,rowDescendants,colDescendants){
  ancs <- unique(rc_relations$ancestor)
  descs <- unique(rc_relations$descendant)
  basal_indexes <- setdiff(ancs,descs) ## ancestors w/o descendant
  
  ### We'll filter rc_table to: (1) only edges with ancestor/descendant relation
  ### i.e. must be in rc_relations
  rc_table <- rc_table[rc_index %in% c(ancs,descs)]
  setkey(rc_table,row.edge,P)
  rc_table <- basal_dominance_filter(rc_tbl=rc_table,rowDescendants)
  return(rc_table)
}

get_lineage <- function(basal_ix,p_thresh=1,rc_relations.=rc_relations,rc_table.=rc_table,row.tree.=row.tree,
                          rowDescendants.=rowDescendants,colDescendants.=colDescendants){
  rct <- rc_table[P<=p_thresh]
  rcm <- rc_relations[max_P<=p_thresh]
  basal_row_edge <- rct[rc_index==basal_ix,row.edge]
  basal_col_edge <- rct[rc_index==basal_ix,col.edge]
  row_descs <- rowDescendants[[basal_row_edge]]
  col_descs <- colDescendants[[basal_col_edge]]
  rc_tbl <- rct[row.edge %in% c(basal_row_edge,row_descs) & 
                       col.edge %in% c(basal_col_edge,col_descs)]
  lineage <- assemble_lineage(rc_tbl,rcm,rt,row.tree)
  lineage[,lineage_id:=basal_ix]
  setkey(lineage,col.edge,P)
  return(lineage)
}


edge_tips <- function(tree){
  grps <- phylofactor::getPhyloGroups(tree)
  mat <- lapply(grps,FUN=function(gg)
    data.table('min'=min(gg[[1]]),'max'=max(gg[[1]])))
  mat <- rbindlist(mat)
  mat[,edge:=1:.N]
  setkey(mat,edge)
  mat[,tip:=(min==max)]
  return(mat)
}

lineage_boxes <- function(lineage,X,colEdgeTips,rowEdgeTips){
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
lstat <- function(id,Lineages,X.=X,colEdgeTips.=colEdgeTips,rowEdgeTips.=rowEdgeTips){
  Xdt <- lineage_boxes(Lineages[lineage_id==id],X,colEdgeTips,rowEdgeTips)
  setkey(Xdt,box)
  Xdt[,box:=factor(box)]
  a <- aov(x~box,data=Xdt)
  ss <-summary(a)[[1]]['box',] %>% unlist
  output <- data.table('lineage_id'=id,'F_stat'=ss['F value'],'msq'=ss['Mean Sq'],'ssq'=ss['Sum Sq'])
  return(output)
}
lineage_stats <- function(Lineages,X,colEdgeTips,rowEdgeTips,cl=NULL){
  lineage_ids <- unique(Lineages$lineage_id)
  fits <- data.table('lineage_id'=lineage_ids,
                     'F_stat'=0)
  if (is.null(cl)){
    for (id in lineage_ids){
      Xdt <- lineage_boxes(Lineages[lineage_id==id],X,colEdgeTips,rowEdgeTips)
      setkey(Xdt,box)
      Xdt[,box:=factor(box)]
      a <- aov(x~box,data=Xdt)
      ss <-summary(a)[[1]]['box',] %>% unlist
      fits[lineage_id==id,F_stat:=ss['F value']]
      fits[lineage_id==id,msq:=ss['Mean Sq']]
      fits[lineage_id==id,ssq:=ss['Sum Sq']]
    }
  } else {
    fits <- parallel::parLapply(cl,lineage_ids,lstat,Lineages=Lineages) %>% rbindlist
  }
  return(fits)
}

global_Fstat <- function(Lineages,X,colEdgeTips,rowEdgeTips){
  lineages <- split(Lineages,f = factor(Lineages$lineage_id)) %>% 
    lapply(lineage_boxes,X,colEdgeTips,rowEdgeTips) %>% rbindlist
  
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

filter_stats <- function(x,Lineages,rc_table,rc_relations,rowDescendants){
  x <- x[order(F_stat,decreasing = T)]
  
  i=1
  while(i<=nrow(x)){
    id <- x$lineage_id[i]
    row.edg <- rc_table[rc_index==id,row.edge]
    descendants <- rowDescendants[[row.edg]]
    ancestors <- edge_ancestors(row.edg,row.tree)
    incompatible_lineages <- setdiff(Lineages[row.edge %in% c(descendants,ancestors),lineage_id],id)
    x <- x[!lineage_id %in% incompatible_lineages]
    i=i+1
  }
  return(x)
}

# edge2node <- function(edges,tree) tree$edge[edges,2]

predict.dendromap <- function(object,...){
  
  lineages <- split(object$Lineages,f = factor(object$Lineages$lineage_id)) %>% 
    lapply(lineage_boxes,object$Data,object$colEdgeTips,object$rowEdgeTips) %>% rbindlist
  
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

get_lineages_and_stats <- function(p_thresh,rc_table,rc_relations,rowDescendants,
                             colEdgeTips,rowEdgeTips,cl=NULL){
  ix_thresh=rc_table[P<=p_thresh,rc_index]
  rct <- rc_table[rc_index %in% ix_thresh]
  rcm <- rc_relations[max_P<=p_thresh]
  basal_indexes <- rcm[ancestor %in% ix_thresh & descendant %in% ix_thresh,setdiff(unique(ancestor),unique(descendant))]
  desc_count <- rcm[ancestor %in% ix_thresh & descendant %in% ix_thresh,list(n=.N),by=ancestor]
  basal_ixs_with_descendants <- intersect(basal_indexes,desc_count[n>1,ancestor])
  if (is.null(cl)){
    Lineages <- lapply(basal_ixs_with_descendants,get_lineage,p_thresh) %>% rbindlist
  } else {
    Lineages <- parallel::parLapply(cl,basal_ixs_with_descendants,get_lineage,p_thresh) %>% rbindlist
  }
  Lineages[,lineage_size:=.N,by=lineage_id]
  Lineages <- Lineages[lineage_size>1]
  stats <- lineage_stats(Lineages,X,colEdgeTips,rowEdgeTips)
  stats <- stats[order(F_stat,decreasing = T)]
  stats <-  filter_stats(stats,Lineages,rct,rcm,rowDescendants)
  Lineages <- Lineages[lineage_id %in% stats$lineage_id]
  
  setkey(Lineages,lineage_id,row.edge)
  setkey(stats,lineage_id)
  
  Lineages <- stats[,c('lineage_id','F_stat')][Lineages]
  return(Lineages)
}

scan_Fstats <- function(pvals,rc_table,rc_relations,rowDescendants,rowEdgeTips,colEdgeTips,cl=NULL){
  ix_prev <- NULL
  Fstats <- rep(NA,length(pvals))
  if (!is.null(cl)){
    ### check cluster?
  }
  start_time <- Sys.time()
  for (i in 1:length(pvals)){
    p_thresh=pvals[i]
    ix_thresh=rc_table[P<=p_thresh,rc_index]
    new_ix <- setdiff(ix_thresh,ix_prev)
    rct <- rc_table[P<=p_thresh]
    rcm <- rc_relations[max_P<=p_thresh]
    
    if (i==1){
      basal_indexes <- rcm[ancestor %in% ix_thresh & descendant %in% ix_thresh,setdiff(unique(ancestor),unique(descendant))]
      desc_count <- rcm[ancestor %in% ix_thresh & descendant %in% ix_thresh,list(n=.N),by=ancestor]
      basal_ixs_with_descendants <- intersect(basal_indexes,desc_count[n>1,ancestor])
      if (is.null(cl)){
        Lineages <- lapply(basal_ixs_with_descendants,get_lineage,p_thresh=p_thresh) %>% rbindlist
      } else {
        Lineages <- parallel::parLapply(cl,basal_ixs_with_descendants,get_lineage,p_thresh=p_thresh) %>% rbindlist
      }
      Lineages[,lineage_size:=.N,by=lineage_id]
      Lineages <- Lineages[lineage_size>1]
    } else {
      ### here, we note that our simplify_rc_table function ensured that any new ix will:
      ### if a descendant row edge from a basal_ix, it cannot have same col.edge as prior basal_ix 
      ##   - this is because P(r,c) < P(r',c) for r' descendant of r <--> basal dominance.
      
      ### Which basal_indexes to we (re)compute?
      ### - ancestors/descendants of new_ix
      ### - note: we may even produce conflicting lineages, e.g. Rheas Gondwanan and Rheas Australian
      ### we DON'T compute: basal_ix NOT on row-tree path (root-tips) from new_ix.
      ### so... 
      ### (1) row_path_ix: Set of ix_thresh on row.tree root-tip path of new_ix
      ### (2) remove all Lineages containing row_path_ix
      ### (3) find basal_ix among row_path_ix
      ### (4) compute Lineages for basal_ix
      ### (5) rbind with old lineages
      ### (6) stats-->filter-->global F stat
      
      prev_row.edges <- rct[rc_index %in% ix_prev,unique(row.edge)]
      new_row.edges <- rct[rc_index %in% new_ix,unique(row.edge)]
      
      ## find new adges ancestral to previous edges
      new_descs <- rowDescendants[new_row.edges]
      new_anc_edges <- sapply(new_descs,FUN=function(a,b) any (b %in% a), b=prev_row.edges) %>% new_row.edges[.]
      
      ## find previous edges ancestral to new edges
      prev_descs <- rowDescendants[prev_row.edges]
      prev_anc_edges <- sapply(prev_descs,FUN=function(a,b) any (b %in% a), b=new_row.edges) %>% prev_row.edges[.]
      
      ## need to recompute all edges which are either equal, ancestral, or descendant to new edges
      row_path_edges <- c(intersect(prev_row.edges,new_row.edges),
                          prev_anc_edges,
                          new_anc_edges) %>% unique
      
      row_path_ix <- rct[row.edge %in% row_path_edges,unique(rc_index)]
      ## NOTE: Since merge_twin_sisters can occasionally introduce a new row.edge not found
      ##       we'll remove all lineage_id's containing row_edges in our row_path_edges set above
      Old_Lineages <- Lineages[!row.edge %in% row_path_edges]
      
      basal_ix <- rcm[ancestor %in% row_path_ix & descendant %in% ix_thresh,setdiff(unique(ancestor),unique(descendant))]
      desc_count <- rcm[ancestor %in% row_path_ix & descendant %in% ix_thresh,list(n=.N),by=ancestor]
      basal_ixs_with_descendants <- intersect(basal_ix,desc_count[n>1,ancestor])
      if (is.null(cl)){
        Lineages <- lapply(basal_ixs_with_descendants,get_lineage,p_thresh=p_thresh) %>% rbindlist
      } else {
        Lineages <- parallel::parLapply(cl,basal_ixs_with_descendants,get_lineage,p_thresh=p_thresh) %>% rbindlist
      }
      Lineages[,lineage_size:=.N,by=lineage_id]
      Lineages <- Lineages[lineage_size>1]
      Lineages <- rbind(Old_Lineages,Lineages)
    }
    
    
    #### need to parallelize lineage_stats
    if (nrow(Lineages)>0){
      stats <- lineage_stats(Lineages,X,colEdgeTips,rowEdgeTips,cl=cl)
      stats <- stats[order(F_stat,decreasing = T)]
      stats <-  filter_stats(stats,Lineages,rct,rcm,rowDescendants)
      Lineages <- Lineages[lineage_id %in% stats$lineage_id]
    }
    if (nrow(Lineages)>0){
      Fstats[i] <- global_Fstat(Lineages,X,colEdgeTips,rowEdgeTips)
    }
    
    ix_prev <- ix_thresh  ## we only have to recompute new_ix that affect our Lineages table
    tm2 <- Sys.time()
    time.elapsed <- signif(difftime(tm2,start_time,units = 'mins'),3)
    GUI.notification <- paste('\r',i,'P-values out of',length(pvals),'scanned in',time.elapsed,'minutes.    ')
    GUI.notification <- paste(GUI.notification,'Estimated time of completion for this step:',
                              as.character(start_time+difftime(tm2,start_time)*length(pvals)/i),
                              '  \r')
    base::cat(GUI.notification)
    utils::flush.console()
  }
  return(data.table('P_thresh'=pvals,'Fstat'=Fstats))
}

dendromap <- function(X,row.tree,col.tree,ncores=NULL,minPval=-Inf,maxPval=0.001,stepsize=10){
  
  # edgeMap, edge_tips for quick descendant calculations -------------------------------------------------------------------
  base::cat('Prepping workspace for edge-based tree traversals')
  rowEdgeMap <- edge_registry(row.tree)
  colEdgeMap <- edge_registry(col.tree)
  row.edges <- 1:Nedge(row.tree)
  col.edges <- 1:Nedge(col.tree)
  rowDescendants <- rowEdgeMap[,list(list(setdiff(seq(edge,edge+n),edge))),by=edge]$V1
  colDescendants <- colEdgeMap[,list(list(setdiff(seq(edge,edge+n),edge))),by=edge]$V1
  colEdgeTips <- edge_tips(col.tree)
  rowEdgeTips <- edge_tips(row.tree)
  
  # edge rc_table -----------------------------------------------------------
  base::cat('\nMaking rc_table')
  rc_table <- make_rc_table(X,row.tree,col.tree,maxPval)
  
  
  # RC_map ------------------------------------------------------------------
  base::cat('\nMaking rc_relations')
  rc_relations <- find_rc_relations(rc_table,rowDescendants,colDescendants)
  
  # simplify rc_table ------------------------------------------------------------------
  base::cat('\nSimplifying rc_table')
  rc_table <- simplify_rc_table(rc_table,rc_relations,rowDescendants,colDescendants)
  rc_relations <- find_rc_relations(rc_table,rowDescendants,colDescendants)
  
  # Scanning F statistics ---------------------------------------------------
  #### set up Pvalues for scanning
  min_pval <-  cbind(rc_table[match(rc_relations$ancestor,rc_index),P],
                     rc_table[match(rc_relations$descendant,rc_index),P]) %>% apply(1,max) %>% min
  rc_relations$max_P <- cbind(rc_table[match(rc_relations$ancestor,rc_index),P],
                       rc_table[match(rc_relations$descendant,rc_index),P]) %>% apply(1,max)
  
  # Prepare cluster -----------------------------------------
  base::cat('\nPreparing cluster')
  if (!is.null(ncores)){
    cl <- parallel::makeCluster(ncores)
    parallel::clusterExport(cl,varlist=c('X','row.tree','col.tree','rc_table','rc_relations',
                                         'rowDescendants','colDescendants','rowEdgeTips','colEdgeTips',
                                         'rowEdgeMap','colEdgeMap'),
                            envir = environment())
    parallel::clusterEvalQ(cl,{library(data.table)
                               library(magrittr)})
    parallel::clusterEvalQ(cl,source('Old_R/edge_dendromap_fcns.R'))
  } else {
    cl <- NULL
  }
  
  pvals <- sort(unique(rc_relations$max_P),decreasing = F) ### need an additional trim
  ps <- rc_table[,list(P=min(P),
                       rc_index=rc_index[which.min(P)]),by=row.edge]
  min_P <- rc_relations[ancestor %in% ps$rc_index,max_P] %>% min
  pvals <- pvals[pvals>=max(min_P,minPval)]
  scanned_pvals <- pvals[seq(1,length(pvals),by=stepsize)]
  base::cat(paste('Beginning coarse scan of',length(scanned_pvals),'P-value thresholds'))
  
  ### Scanning
  Fscan <- tryCatch(scan_Fstats(scanned_pvals,rc_table,rc_relations,rowDescendants,
                               rowEdgeTips,colEdgeTips,cl=cl),
                    error=function(e) e)
  if (!'data.table' %in% class(Fscan)){
    if (!is.null(cl)){
      parallel::stopCluster(cl)
      rm('cl')
      gc()
    }
    stop(Fscan)
  }
  
  
  ### Refining - scan every pval around maximum
  ix <- which(scanned_pvals==Fscan[Fstat==max(Fstat),P_thresh])
  if (ix<=2){
    refined_pvals <- pvals[pvals<=scanned_pvals[4]]
  } else if (ix>=(length(scanned_pvals)-2)){
    refined_pvals <- pvals[pvals>=scanned_pvals[length(scanned_pvals)-3]]
  } else {
    refined_pvals  <- pvals[pvals>=scanned_pvals[ix-2] & pvals<=scanned_pvals[ix+2]]
  }
  base::cat(paste('Beginning refined scan of',length(refined_pvals),'P-value thresholds'))
  
  Fscan_refined <- tryCatch(scan_Fstats(refined_pvals,rc_table,rc_relations,rowDescendants,
                                       rowEdgeTips,colEdgeTips,cl=cl),
                            error=function(e) e)
  if (!'data.table' %in% class(Fscan_refined)){
    if (!is.null(cl)){
      parallel::stopCluster(cl)
      rm('cl')
      gc()
    }
    stop(Fscan_refined)
  }
  Fscan <- rbind(Fscan,Fscan_refined)
  Fscan <- Fscan[!duplicated(Fscan)]
  
  # Preparing Output ------------------------------------------------------
  
  Lineages <- Fscan[Fstat==max(Fstat),P_thresh] %>% 
    get_lineages_and_stats(rc_table,rc_relations,rowDescendants,
                     colEdgeTips,rowEdgeTips,cl)
  
  
  if (!is.null(cl)){
    parallel::stopCluster(cl)
    rm('cl')
    gc()
  }
  
  setkey(Lineages,row.edge)
  edge2node <- function(edge,tree) tree$edge[edges,2]
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
                 'rowDescendants'=rowDescendants,
                 'colDescendants'=colDescendants,
                 'F_scan'=Fscan)
  class(object) <- 'dendromap'
  return(object)
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
  desc_edges <- dm$rowDescendants[[basal_edge]]
  tips <- phangorn::Descendants(dm$row.tree,dm$row.tree$edge[basal_edge,2],'tips')[[1]]
  rt <- ape::drop.tip(dm$row.tree,setdiff(dm$row.tree$tip.label,dm$row.tree$tip.label[tips]))
  rownames(rt$edge) <- desc_edges
  
  L[row.edge!=basal_edge,row.node:=rt$edge[as.character(row.edge),2]]
  L[row.edge==basal_edge,row.node:=(length(rt$tip.label)+1)]
  
  temp_dm <- list('Data'=dm$Data[row.tree$tip.label[tips],],
                  'Lineages'=L,
                  'col.tree'=dm$col.tree,
                  'row.tree'=rt)
  class(temp_dm) <- 'dendromap'
  pl <- plot.dendromap(temp_dm,...)
  return(pl)
}
