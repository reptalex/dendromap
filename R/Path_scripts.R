##Here I give three functions that return:
## 1) all the non-terminal paths starting from the root
## 2) all the depth intervals of these paths
## 3) the orientation of all those paths
## The outputs of all three functions are lists, and the element numbers in these lists correspond to the same path 
## The depths can be used to calculate Gillepsie propensities and manifest events based on these propensities
## The depths can also be used to calculate null distributions for paths of interest, to create confidence statstics that check for false positive patterns between the row and the column tree

library(ape)
library(phangorn)

trimpaths <- function(Tree) {
  Nodepath <- nodepath(Tree, from=Ntip(tree)+1)
  trimedpath <- list(length=length(Nodepath))
  for (i in 1:(length(Nodepath))) {
    trimedpath[[i]] <- Nodepath[[i]][-length(Nodepath[[i]])]
  }
  trimedpath <- unique(trimedpath)
  return(trimedpath)
}

pathdepths <- function(trimedpath, tree) {
  Depths <- node.depth.edgelength(tree)
  for (j in 1:(length(trimedpath))) {
    for (i in 1:length(trimedpath[[j]])) {
      for (k in (Ntip(tree)+1):length(Depths)) {
        trimedpath[[j]] <- replace(trimedpath[[j]], trimedpath[[j]]==k, Depths[k])
      }
    }
  }
  return(trimedpath)
}

pathorient <- function(trimedpath, tree) {
  Desc <- Descendants(tree, (Ntip(tree)+1), type="all")
  desc <- vector(length = length(Desc))
  is.even <- function(x) x %% 2 == 0
  for (i in 1:(length(Desc))) {
    if (isTRUE(is.even(i))) {
      desc[i] <- 'a'
    } else {
      desc[i] <- 'b'
    }
  }
  
  for (j in 1:(length(trimedpath))) {
    for (i in 1:length(trimedpath[[j]])) {
      for (k in (Ntip(tree)+1):length(desc)) {
        trimedpath[[j]] <- replace(trimedpath[[j]], trimedpath[[j]]==k, desc[k])
      }
    }
  }
  return(trimedpath)
}
