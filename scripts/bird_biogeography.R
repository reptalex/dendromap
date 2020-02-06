library(dendromap)


load('data/birds/bird_conts.Rdata')
load('data/birds/PTrees.Rdata')
row.tree <- PTrees[[1]]
col.tree <- readRDS('data/birds/continental_tree')

colnames(bird_conts)[c(2,5,6)] <- c('Eurasia','NAmerica','SAmerica')

# checking ----------------------------------------------------------------

X <- bird_conts[,2:6] %>% as.matrix
rownames(X) <- bird_conts$species
colnames(X) <- colnames(bird_conts)[2:6]

all(rownames(X) %in% row.tree$tip.label)
all(colnames(X) %in% col.tree$tip.label)


row.tree <- drop.tip(row.tree,setdiff(row.tree$tip.label,rownames(X)))

X <- X[row.tree$tip.label,col.tree$tip.label]


dm <- dendromap2(X,row.tree,col.tree,ncores=7,estimate_runtime = TRUE)
