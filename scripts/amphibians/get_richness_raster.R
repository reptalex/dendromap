# This script defines a function to extract a raster of species richness for an arbitrary list of amphibian species

# The function takes two inputs: a list of species names and the filepath to the locations of the single-species rasters (e.g.
# as output by make_rasters.R)

get_richness_raster <- function(species_list, raster_location){
  if(length(species_list) < 1){stop('species_list must have at least one element')}
  if(length(species_list) != length(unique(species_list))){warning('species_list contains duplicate elements')}
  file_list <- list.files(raster_location, full.names = T)
  grd_files <- file_list[grepl(".grd", file_list)]
  
  raster_i <- raster::raster(grd_files[[1]])
  richness_raster <- raster_i
  
  if(length(species_list) > 1){
    for(i in 2:length(species_list)){
      raster_i <- raster::raster(grd_files[[i]])
      richness_raster <- richness_raster + raster_i
    }
  }
  return(richness_raster)
}