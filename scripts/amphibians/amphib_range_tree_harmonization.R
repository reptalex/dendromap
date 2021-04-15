# This script:
# 1. loads the phylogeny from https://datadryad.org/stash/dataset/doi:10.5061/dryad.cc3n6j5
# 2. the amphibian maps from IUCN/Natureserve
# 3. loads a shapefile of continents including major land-bridge islands
# 4. subsets the range maps to those that intersect the continents
# 5. where possible, automatically standardizes names from the phylogeny and the range maps to GBIF taxonomy
#         NOTE: THIS INTRODUCES SOME DUPLICATION IN THE TIP LABELS OF THE PHYLOGENY WHICH NEEDS TO BE RESOLVED SOMEHOW
# 6. manually for all range-maps labels that could not be standardized, or are otherwise absent from the 
#         standardized phylogeny tip labels, manually harmonizes taxonomy or drops the species from the range
#         map dataset as appropriate
# 7. saves the reduced, harmonized range data (RDS save of sf object) and the harmonized phylogeny (RDS save of phylo object)

library(sf)
library(ape)
library(ggplot2)
`%ni%` <- Negate(`%in%`)

# Amphibian phylogeny
amphib_phylo <- read.tree('/Users/jacobsocolar/Downloads/doi_10.5061_dryad.cc3n6j5__v1/amph_shl_new_Consensus_7238.tre')
scinames <- amphib_phylo$tip.label

# Amphibian ranges
amphib_ranges1 <- read_sf('/Users/jacobsocolar/Downloads/AMPHIBIANS/AMPHIBIANS.shp')

# Continents (restrict to mainlands plus large land-bridge islands)
setwd("/Users/JacobSocolar/Dropbox/Work/Diversity_accum")
load("continents.Rdata")
geometry <- lapply(continents, function(x) st_as_sf(x))
conts <- do.call(rbind, geometry)
conts$id <- names(continents)
conts2 <- st_transform(conts, st_crs(amphib_ranges1))
ggplot(conts2) + geom_sf()

# Which amphibians intersect with continental mainlands-plus-land-bridge-islands
amphib_conts <- as.data.frame(st_intersects(amphib_ranges1, conts2, sparse = F))
names(amphib_conts) <- names(continents)

# Reduce amphibian ranges to continental species
amphib_ranges <- amphib_ranges1[rowSums(amphib_conts[,2:6]) > 0, ]

# Standardize the species names from range data to GBIF
range_standardized1 <- traitdataform::get_gbif_taxonomy(amphib_ranges$binomial[1:2552], subspecies = F, verbose = T)
range_standardized2 <- traitdataform::get_gbif_taxonomy(amphib_ranges$binomial[2554:nrow(amphib_ranges)], subspecies = F, verbose = T)
range_standardized <- c(range_standardized1$scientificName, NA, range_standardized2$scientificName)
range_extras <- amphib_ranges$binomial[is.na(range_standardized)]

# Standardize the species names from phylogeny to GBIF
tree_standardized1 <- traitdataform::get_gbif_taxonomy(scinames[1:5287], subspecies = F, verbose = T)
tree_standardized2 <- traitdataform::get_gbif_taxonomy(scinames[5290:6110], subspecies = F, verbose = T)
tree_standardized3 <- traitdataform::get_gbif_taxonomy(scinames[6112:length(scinames)], subspecies = F, verbose = T)
tree_standardized <- c(tree_standardized1$scientificName, rep(NA, 2), tree_standardized2$scientificName, NA, tree_standardized3$scientificName)
tree_extras <- scinames[is.na(tree_standardized)]

# Reinclude names that could not be standardized
tree_standardized[is.na(tree_standardized)] <- scinames[is.na(tree_standardized)]
range_standardized[is.na(range_standardized)] <- amphib_ranges$binomial[is.na(range_standardized)]

# Check species lists against one another
tree_standardized[tree_standardized %ni% range_standardized] # Most of these species are on islands

range_standardized[range_standardized %ni% tree_standardized]
# This is more species that didn't get standardized than I'd hope.  Check how many of them are 
# recently described:
range_authors <- c(range_standardized1$author, NA, range_standardized2$author)
authors2 <- range_authors[range_standardized %ni% tree_standardized]
range_dates <- as.numeric(stringr::str_extract(authors2, "[0-9]{4}"))
hist(range_dates)
sum(range_dates > 2013, na.rm = T)
sum(range_dates <= 2013, na.rm = T)
# The overwhelming majority are 2014 or later.  We'll assume that these are likely to be omitted
# from the phylogeny because they're too recently described.  Therefore, we omit from analysis.
new_species <- as.numeric(stringr::str_extract(range_authors, "[0-9]{4}")) >= 2014
include_species <- (!new_species) | (range_standardized %in% tree_standardized)
include_species[is.na(include_species)] <- T


# Now we build a version of amphib_ranges with names standardized to the tree:
amphib_ranges$name_standardized <- range_standardized

# Get rid of the newly described species
amphib_ranges_reduced <- amphib_ranges[include_species,]

# Update all remaining conflicts by hand
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Nanohyla marmorata"] <- "Microhyla marmorata"
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Feihyla samkosensis",]  # Newly described
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Afrixalus crotalus"] <- "Afrixalus aureus" # Split
tree_standardized[scinames == "Leptopelis_barbouri"] <- "Leptopelis grandiceps" # Weird taxonomic shuffling
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Oreolalax sterlingae",]  # Newly described
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Rohanixalus vittatus"] <- "Feihyla vittata"
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Eleutherodactylus campi"] <- "Eleutherodactylus cystignathoides"
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Telmatobufo ignotus",]  # Apparently unmapped
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Cardioglossa inornata",]  # Newly described
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Pachytriton xanthospilos",]  # Newly-ish described
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Chiropterotriton perotensis",]  # Newly described
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Nanohyla annectens"] <- "Microhyla annectens"
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Rohanixalus punctatus"] <- "Chiromantis punctatus"
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Hyla imitator"] <- "Dendropsophus_imitator"
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Lithobates johni",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Cynops yunnanensis",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Taricha sierra"] <- "Taricha sierrae"
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Hynobius vandenburghi",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Zhangixalus leucofasciatus",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Nanohyla arboricola",]  # Newly described
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Pristimantis imthurni",]  # Newly described
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Hyperolius papyri",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Nanohyla nanapollexa"] <- "Microhyla nanapollexa"
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Thorius maxillabrochus",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Gastrotheca lojana",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Indirana tenuilingua"] <- "Indirana_tenuilingua"
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Ololygon skuki",]  # Newly described
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Theloderma albopunctatum",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Theloderma baibungense",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Callulops wondiwoiensis",]  # Newly described
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Hyloxalus peruvianus"] <- "Allobates_peruvianus"
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Humerana oatesii"] <- "Hylarana_oatesii"
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Xenopus calcaratus",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Lithobates palustris"] <- "Rana_palustris"
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Bijurana nicobariensis"] <- "Amnirana nicobariensis"
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Chiropterotriton totonacus",]  # Newly described
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Amolops chayuensis",]  # Newly described
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Nanorana sichuanensis",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Gracixalus quyeti",]  # Newly-ish described
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Oreophryne cameroni",]  # Newly described
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Rohanixalus baladika",]  # Newly described
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Nanohyla pulchella",]  # Newly described
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Kalophrynus barioensis",]  # Newly described
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Chiropterotriton casasi",]  # Newly described
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Nanohyla annamensis"] <- "Microhyla annamensis"
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Proceratophrys salvatori"] <- "Odontophrynus salvatori"
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Isthmura sierraoccidentalis",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Rohanixalus hansenae"] <- "Feihyla hansenae"
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Craugastor blairi",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Nanohyla perparva"] <- "Microhyla perparva"
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Tepuihyla warreni"] <- "Dendropsophus_warreni"
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Afrixalus brachycnemis",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Humerana humeralis"] <- "Hylarana_humeralis"
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Boana claresignata"] <- "Bokermannohyla claresignata"
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Rohanixalus nauli",]  # Newly described
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Rhinella beebei",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Anomaloglossus megacephalus",]  # Newly described
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Eleutherodactylus orarius",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Xenopus poweri",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Gracixalus carinensis"] <- "Gracilixalus_carinensis"
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Nyctimystes tyleri"] <- "Nyctimystes_michaeltyleri"
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Nanohyla petrigena"] <- "Microhyla petrigena"
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Rohanixalus shyamrupus"] <- "Chiromantis shyamrupus"
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Oreophryne parkopanorum",]  # Newly described
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Pristimantis jamescameroni",]  # Newly described
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Trachycephalus typhonius"] <- "Trachycephalus venulosus"
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Chiropterotriton melipona",]  # Newly described
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Eleutherodactylus bilineatus"] <- "Noblella_bilineata"
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Humerana miopus"] <- "Hylarana_miopus"
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Rohanixalus marginis",]  # Newly described
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Petropedetes newtonii",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Callulops neuhaussi",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Bolitoglossa la",]  # Newly described
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Chiropterotriton ceronorum",]  # Newly described
amphib_ranges_reduced$name_standardized[amphib_ranges_reduced$name_standardized == "Boana clepsydra"] <- "Bokermannohyla clepsydra"
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Austrochaperina punctata",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Dendropsophus baileyi",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Pristimantis ameliae",]  # Newly described
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Dendropsophus baileyi",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Sclerophrys pusilla",]  # Apparent lump; closely related species all share same continental association
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Ololygon atrata",]  # Missing; unexplained
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Ameerega boehmei",]  # Missing; unexplained
amphib_ranges_reduced <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized != "Ptychadena boettgeri",]  # Missing; unexplained

saveRDS(amphib_ranges_reduced, "/Users/jacobsocolar/Dropbox/Work/amphib_phylo_richness/amphib_ranges_reduced.RDS")


amphib_phylo2 <- amphib_phylo
amphib_phylo2$tip.label <- tree_standardized
saveRDS(amphib_phylo2, "/Users/jacobsocolar/Dropbox/Work/amphib_phylo_richness/amphib_phylogeny2.RDS")
