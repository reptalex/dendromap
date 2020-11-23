## load dendromap
# library(biomformat)
# library(rbiom)
library(ape)
library(data.table)
library(magrittr)
row.tree <- read.tree('data/tetrapod_microbiome_analysis-master/tree.nwk')
col.tree <- read.tree('data/tetrapod_microbiome_analysis-master/trees/total_timetree_names.aves.nwk.tre')

# 
dat <- rbiom::read.biom('data/tetrapod_microbiome_analysis-master/aves.biom',tree=FALSE)


# filtering ---------------------------------------------------------------
seqs <- rownames(dat$taxonomy)
all(seqs %in% row.tree$tip.label)
row.tree <- drop.tip(row.tree,setdiff(row.tree$tip.label,seqs))

md <- read.csv('data/tetrapod_microbiome_analysis-master/metadata/5per.2.18.19.short.txt',sep='\t') %>% as.data.table
md[,species:=gsub(' ','_',Species_name)]

col.tree <- drop.tip(col.tree,setdiff(col.tree$tip.label,md$species))
md <- md[match(col.tree$tip.label,species),]
all.equal(md$species,col.tree$tip.label)

X <- rbiom::counts(dat)

birds <- md[SampleID %in% colnames(X),species]

col.tree <- drop.tip(col.tree,setdiff(col.tree$tip.label,birds))


X <- X[row.tree$tip.label,md[species %in% birds,SampleID]]

zeros <- rownames(X)[which(rowSums(X)==0)]

row.tree <- drop.tip(row.tree,zeros)
X <- X[row.tree$tip.label,]

any(colSums(X) != 1e4)

mins <- apply(X,1,FUN=function(y) min(y[y>0]))
## we'll keep sequences with over 1/1e3 prevalence, i.e. minimum 10 counts.

keepers <- names(mins)[mins>=10]

row.tree <- drop.tip(row.tree,setdiff(row.tree$tip.label,keepers))
X <- X[row.tree$tip.label,]

colnames(X) <- col.tree$tip.label
dm <- dendromap(X,row.tree,col.tree,ncores=7)
save(list=ls(),file='data/tetrapod_microbiome_analysis-master/dendromap_workspace')



tinamou <-col.tree$tip.label[grepl('Nothura',col.tree$tip.label)]
rhea <- col.tree$tip.label[grepl('Rhea',col.tree$tip.label)]
anc <- ape::getMRCA(col.tree,which(col.tree$tip.label %in% c(rhea,tinamou)))


dm$Lineages[,basal_rc:=min(col.node),lineage_id]
dm$Lineages[col.node==anc] ## 6 Lineages contain the Rheas! all are basal_rc
rhea_containing_lineages <- dm$Lineages[col.node==387,Lineage]
dm$Lineages[Lineage %in% rhea_containing_lineages]


L_sum <- dm$Lineages[,list(n_statements=.N,
                           basal_row=min(row.node),
                           basal_col=min(col.node)),by=Lineage][order(n_statements)]

### e.g. 6326

lineage_plot(dm,id=6326)
plot.dendromap(dm)
