# rm(list=ls())
library(ape)
library(phytools)
library(data.table)
library(maptools)
library(raster)

# Making map for introduction

#########################
#all.life.forms <- function(reference_table, all_vars, scoring, chunks=10000) {
#  result_traits <- data.frame(species=reference_table$wcvp_name, life_form=NA)
#  log_preliminary <- seq(1, nrow(reference_table), chunks)[-1]
#  for(species_index in 1:nrow(reference_table)) {
#    if(species_index %in% log_preliminary) {
#      save(result_traits, file="result_traits_tmp.Rsave")
#    }
#    cat(species_index, "\r")
#    wcvp_subset <- subset(all_vars, all_vars$taxon_name == reference_table$wcvp_name[species_index])
#    life_form <- tail(names(sort(table(wcvp_subset$lifeform_description))),1)
#    if(life_form=="") {
#      # try second one
#      life_form <- tail(names(sort(table(wcvp_subset$lifeform_description))),2)
#      if(life_form[1]=="") {
#        result_traits[species_index,2] <- "no_life_form_on_database"
#      } else {
#        result_traits[species_index,2] <- scoring$scoring1[which(scoring$all_life_forms==life_form[1])]
#      }
#      } else { 
#      result_traits[species_index,2]  <- scoring$scoring1[which(scoring$all_life_forms==life_form)]
#      }
#    }  
#  return(result_traits)
#}
#########################
organize.bubble.plot <- function(trait_table, reference_table, all_vars, twgd_data) {
  tmp_reference_table <- subset(reference_table, reference_table$wcvp_name %in% unique(trait_table$species))
  wcvp_subset <- subset(all_vars, all_vars$taxon_name %in% tmp_reference_table$wcvp_name)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$introduced==0)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$extinct==0)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$location_doubtful==0)
  
  focal_areas <- unique(wcvp_subset$area_code_l3)
  results <- matrix(nrow=0, ncol=5)
  for(i in 1:length(focal_areas)) {
    one_area <- focal_areas[i]
    one_subset <- subset(wcvp_subset, wcvp_subset$area_code_l3==one_area)
    sp_rich <- length(unique(one_subset$taxon_name))
    family_rich <- length(unique(one_subset$family))
    area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% one_area),]
    if(nrow(area_plus_buffer)>0) {
      centroids <- rgeos::gCentroid(area_plus_buffer, byid=TRUE)
      lon <- extent(centroids)[1]
      lat <- extent(centroids)[3]
      results <- rbind(results, cbind(sp_rich, family_rich, one_area, lon, lat))
    }
    cat(i, "\r")
  }
  results <- as.data.frame(results)
  results$sp_rich <- as.numeric(results$sp_rich)
  results$family_rich <- as.numeric(results$family_rich)
  results$lon <- as.numeric(results$lon)
  results$lat <- as.numeric(results$lat)
  return(results)
}
#########################
sum.twgd.trait <- function(one_dataset,twgd_data,all_vars) {
  tmp_rasters <- list()
  for(i in 1:nrow(one_dataset)) {
    wcvp_subset <- subset(all_vars, all_vars$taxon_name == one_dataset$species[i])
    occ_areas <- wcvp_subset$area_code_l3
    area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% occ_areas),]
    tmp0 <- raster(area_plus_buffer, res=0.1)
    tmp0[]  <- 1  
    tmp1 <- crop(tmp0, extent(area_plus_buffer))
    tmp2 <- mask(tmp1, area_plus_buffer)
    tmp_rasters[[i]] <- tmp2
    cat(paste0(i, " out of ", nrow(one_dataset)), "\r")
  }
  return(tmp_rasters)
}

#########################
#transform.raster <- function(organized_table_for_plot1, twgd_data) {
#  template.map <- readRDS("figures/global_map_life_form/template.map.Rdata")
#  tmp.raster.list <- list()
#  final_map <- template.map
#  for(i in 1:nrow(organized_table_for_plot1)) {
#    one_area <- organized_table_for_plot1$one_area[i]
#    area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) == one_area),]
#    if(nrow(area_plus_buffer)>0) {
#      final_map <- rasterize(area_plus_buffer, final_map, field = organized_table_for_plot1$sp_rich[i], update = TRUE)
#      cat(i, "\r")
#    }
#  }
#  return(final_map)  
#}
#########################
# inefficient:
#transform.raster <- function(organized_table_for_plot1, twgd_data) {
#  template.map <- readRDS("figures/global_map_life_form/template.map.Rdata")
#  tmp.raster.list <- list()
#  for(i in 1:nrow(organized_table_for_plot1)) {
#    one_area <- organized_table_for_plot1$one_area[i]
#    area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) == one_area),]
#    if(nrow(area_plus_buffer)>0) {
#      tmp0 <- raster(area_plus_buffer, res=0.1)
#      tmp0[]  <- organized_table_for_plot1$sp_rich[i]
#      tmp1 <- crop(tmp0, extent(area_plus_buffer))
#      tmp2 <- mask(tmp1, area_plus_buffer)
#      tmp.raster.list[[i]] <- tmp2
#      cat(i, "\r")
#    }
#  }
#  cat("\n")
#  cat("calculating species richness...")
#  tmp.raster.list[which(unlist(lapply(tmp.raster.list, function(x) extent(x)[3] < -90)))] <- NULL
#  tmp.raster.list[which(unlist(lapply(tmp.raster.list, function(x) extent(x)[4] < -60)))] <- NULL
#  tmp.raster.list0 <- lapply(tmp.raster.list, resample, template.map)
#  tmp.raster.list0[which(unlist(lapply(tmp.raster.list0, function(x) all(is.na(x[])))))] <- NULL
#  for(u in 1:length(tmp.raster.list0)) {
#    r0 <- tmp.raster.list0[[u]]
#    r0[which(is.na(r0[]))] <- 0
#    tmp.raster.list0[[u]] <- r0
#  }
#  sprichness_map <- raster::calc(raster::stack(tmp.raster.list0), sum)
#  return(sprichness_map)  
#}
#########################
#########################
# If local
#setwd("~/Desktop/WCVP_special_issue/James_perennial_annual/life_history_houwie")
source("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/WCVPtools_functions.R")
dist_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
#-----------------------------
# If labcomputer
# setwd("~/life_history_houwie")
#source("../WCVPtools/WCVPtools_functions.R")
#dist_sample <- read.table("/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
#names_sample <- read.table("/wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

#########################
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")

#########################
# reference table for taxized names
#-----------------------------
# If local
reference_table <- list.files("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/taxized_reference_tables", full.names = T)
reference_table <- do.call(rbind, lapply(reference_table, read.csv))
#-----------------------------
# If labcomputer
#reference_table <- list.files("../WCVPtools/taxized_reference_tables", full.names = T)
#reference_table <- do.call(rbind, lapply(reference_table, read.csv))
#-----------------------------
#scoring = read.csv("2022-03-19_life_form.csv")

#########################
#life_forms <- all.life.forms(reference_table[1:100,], all_vars, scoring)
#life_forms <- subset(life_forms, life_forms$life_form!="no_life_form_on_database")
#write.csv(life_forms, file="figures/global_map_life_form/all_life_forms.csv", row.names=F)

###########
# Looking at the WCVP table and TDWG to clean GBIF points
#-----------------------------
# If local
path="/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/wgsrpd-master/level3/level3.shp"
#-----------------------------
# If labcomputer
# path="../WCVPtools/wgsrpd-master/level3/level3.shp"
#-----------------------------

twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))
#life_forms <- read.csv("figures/global_map_life_form/all_life_forms_cleaned.csv")
#organized_table_for_plot1 <- organize.bubble.plot(life_forms, reference_table, all_vars, twgd_data)


#total_richness <- sprichness_map
#writeRaster(total_richness, file="total_richness.tif")

life_forms <- read.csv("figures/global_map_life_form/all_life_forms_cleaned.csv")
organized_table_for_plot_total <- organize.bubble.plot(life_forms, reference_table, all_vars, twgd_data)
organized_table_for_plot_annuals <- organize.bubble.plot(subset(life_forms, life_forms$life_form=="annual"), reference_table, all_vars, twgd_data)
#organized_table_for_plot_perennials <- organize.bubble.plot(subset(life_forms, life_forms$life_form=="perennial"), reference_table, all_vars, twgd_data)

nas <- rep(NA, nrow(organized_table_for_plot_annuals))
proportion_table <- data.frame(sp_rich_prop=nas, sp_rich_total=nas, one_area=nas, lon=nas, lat=nas)
for(i in 1:nrow(organized_table_for_plot_annuals)) {
  annual_sp_rich <- organized_table_for_plot_annuals$sp_rich[i]
  one_area <- organized_table_for_plot_annuals$one_area[i]
  total_sp_rich <- organized_table_for_plot_total$sp_rich[organized_table_for_plot_total$one_area == one_area]
  one_proportion <- round(annual_sp_rich / total_sp_rich, 3)
  proportion_table$sp_rich_prop[i] <- one_proportion
  proportion_table$sp_rich_total[i] <- total_sp_rich
  proportion_table$one_area[i] <- one_area
  proportion_table$lon[i] <- organized_table_for_plot_annuals$lon[i]
  proportion_table$lat[i] <- organized_table_for_plot_annuals$lat[i]
}

library(ggplot2)
library(maps)
library(ggthemes)
library(viridis)

twgd_data01 <- sf::st_as_sf(twgd_data)
twgd_data01 <- merge(twgd_data01, proportion_table, by.x="LEVEL3_COD", by.y="one_area")

# portrait, 9 x 7
tmp_map1 <- ggplot(data = twgd_data01) +
  geom_sf(aes(fill = sp_rich_prop, lwd=0.5)) +
  scale_fill_viridis_c(option = "plasma", trans="sqrt") +
  theme_classic() 


tmp_map2 <- ggplot(data = twgd_data01) +
  geom_sf(aes(fill = sp_rich_total)) +
  scale_fill_viridis_c(option = "plasma", trans="sqrt") +
  theme_classic()

mybreaks = seq(0,0.6, by=0.1)
tmp_map3 <- tmp_map2 +
  geom_point(data = twgd_data01, aes(x=lon, y=lat, size=sp_rich_prop), alpha=0.5, shape=20, stroke=FALSE) +
  scale_size_continuous(name="proportion annuals", range=c(0,20), breaks=c(0.05,0.1, 0.5)) 
#+
  #scale_color_viridis(option="plasma", alpha=0.25) 

# alpha and color for diversity and proportion
# using raster:
proportion_raster <- transform.raster(proportion_table, twgd_data)
