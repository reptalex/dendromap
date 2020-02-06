library(phylofactor)
tr <- NULL
tr$tip.label <- c('Eurasia','NAmerica','Australia','Africa','SAmerica')
tr$edge <- c(6,7,
             7,1,
             7,2,
             6,8,
             8,3,
             8,9,
             9,4,
             9,5) %>%
              matrix(byrow=T,ncol=2)
tr$edge.length <- rep(1,nrow(tr$edge))
tr$Nnode <- 4
class(tr) <- 'phylo'

saveRDS(tr,'data/birds/continental_tree')