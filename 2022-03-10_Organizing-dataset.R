##
# setwd("~/Desktop/James_perennial_annual/life_history_houwie")
# rm(list=ls())
library(ape)
library(data.table)
library(rgbif)
library(maptools)
library(raster)
#data("wrld_simpl")

FilterWCVP <- function(points, all_vars, reference_table, species= "scientificName", lon="decimalLongitude", lat="decimalLatitude", path="wgsrpd-master/level3/level3.shp") {
  npoints_start <- nrow(points)
  tmp_points = as.data.frame(points)
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  # Load shape files and make sure they have the same name as the WCVP column with the TDWG areas
  twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))
  dubiousGBIF_ids <- c()
  for(species_index in 1:nrow(reference_table)) {
    gbif_subset <- subset(tmp_points, tmp_points$scientificName == reference_table$gbif_name[species_index])
    if(nrow(gbif_subset)!=0) {
      wcvp_subset <- subset(all_vars, all_vars$taxon_name == reference_table$wcvp_name[species_index])
      occ_areas <- wcvp_subset$area_code_l3
      area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% occ_areas),]
      coords <- gbif_subset[,c("x","y")]
      sp::coordinates(coords) <- ~ x + y
      answer <- which(is.na(sp::over(coords, area_plus_buffer)[,3]))
      if(length(answer) != 0) {
        dubiousGBIF_ids <- c(dubiousGBIF_ids, as.character(gbif_subset$gbifID[answer]))
      }
    }
  }
  cleaned_points <- subset(points, !as.character(points$gbifID) %in% dubiousGBIF_ids)
  npoints_end <- nrow(cleaned_points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(cleaned_points)
}

# Getting climate means per TWGD region
MeanRasterWCVP <- function(path_raster="3_Landscape_instability/bio_1_instability.tif", path_tdwg="wgsrpd-master/level3/level3.shp") {
  twgd_data <- suppressWarnings(maptools::readShapeSpatial(path_tdwg))
  raster_example <- raster(path_raster)
  template_map <- raster_example
  template_map[!is.na(template_map[])] <- 0
  template_map <- aggregate(template_map, fact=25)
  raster_list <- list()
  for(area_index in 1:length(twgd_data)){
    one_area <- twgd_data[area_index,]
    cropped_raster <- NULL
    try(cropped_raster <- mask(crop(raster_example, one_area), one_area))
    if(!is.null(cropped_raster)) {
      mean_value <- mean(subset(cropped_raster[], !is.na(cropped_raster[])))
      template_raster <- cropped_raster
      template_raster[!is.na(template_raster[])] <- mean_value
      template_raster <- aggregate(template_raster, fact=25)
      template_raster <- raster::resample(template_raster, template_map)
      raster_list[[area_index]] <- raster::mask(template_raster, template_map) 
      print(area_index)
    } else { next }
  }
  raster_list[which(unlist(lapply(raster_list, is.null)))] <- NULL
  mm <- do.call(merge, raster_list)
  return(mm)
}



gbif_data <- fread("0106966-210914110416597.csv") # load the table you downloaded from GBIF

# Looking at the WCVP table and TDWG to clean GBIF points
cleaned_gbif <- FilterWCVP(gbif_data, all_vars, reference_table) # This will filter the GBIF points acording to WCVP


#------------------------------------ still developing this:
test <- MeanRasterWCVP(path_raster="3_Landscape_instability/bio_1_instability.tif", path_tdwg="wgsrpd-master/level3/level3.shp")
