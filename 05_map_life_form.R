# send names to gbif
setwd("~/Desktop/WCVP_special_issue/James_perennial_annual/life_history_houwie")
# rm(list=ls())
library(ape)
library(phytools)
library(data.table)
library(maptools)
library(raster)

#########################
all.life.forms <- function(reference_table, all_vars, scoring) {
  result_traits <- data.frame(species=reference_table$wcvp_name, life_form=NA)
  for(species_index in 1:nrow(reference_table)) {
    cat(species_index, "\r")
    wcvp_subset <- subset(all_vars, all_vars$taxon_name == reference_table$wcvp_name[species_index])
    life_form <- tail(names(sort(table(wcvp_subset$lifeform_description))),1)
    if(life_form=="") {
      # try second one
      life_form <- tail(names(sort(table(wcvp_subset$lifeform_description))),2)
      if(life_form[1]=="") {
        result_traits[species_index,2] <- "no_life_form_on_database"
      } else {
        result_traits[species_index,2] <- life_form[1] 
      }
      } else { 
      result_traits[species_index,2]  <- life_form
      }
    }  
  for(j in sequence(nrow(result_traits))) {
    life_form2 <- result_traits[j,2]
    if(life_form2!="no_life_form_on_database") {
      result_traits[j,2] <- scoring$scoring1[which(scoring$all_life_forms==life_form2)]        
    }
  } 
  return(result_traits)
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
GetSpRichness <- function (ranges) {
  template.map <- readRDS("global_map_life_form/template.map.Rdata")
  tmp.raster.list <- list()
  for (i in 1:length(ranges)) {
    r1 <- ranges[[i]]
    r1 <- raster::resample(r1, template.map)
    r1[is.na(r1)] <- 0
    tmp.raster.list[[i]] <- raster::mask(r1, template.map)
    cat(paste0(i, " out of ", length(ranges)), "\r")
  }
  cat("\n")
  names(tmp.raster.list) <- names(ranges)
  cat("calculating species richness...")
  sprichness_map <- raster::calc(raster::stack(tmp.raster.list), sum)
  return(sprichness_map)
}
#########################
#########################
#########################
dist_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")

#########################
reference_table <- list.files("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/taxized_reference_tables", full.names = T)
reference_table <- do.call(rbind, lapply(reference_table, read.csv))
scoring = read.csv("2022-03-19_life_form.csv")

#########################
life_forms <- all.life.forms(reference_table[1:100,], all_vars, scoring)
life_forms <- subset(life_forms, life_forms$life_form!="no_life_form_on_database")
write.csv(life_forms, file="global_map_life_form/all_life_forms.csv", row.names=F)

#########################
path_tdwg = "/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools//wgsrpd-master/level3/level3.shp"
twgd_data <- suppressWarnings(maptools::readShapeSpatial(path_tdwg))
annuals <- subset(life_forms, life_forms$life_form=="annual")
perennials <- subset(life_forms, life_forms$life_form=="perennial")

#########################
annuals_raster <- sum.twgd.trait(annuals,twgd_data,all_vars)
save(annuals_raster, file="global_map_life_form/annuals_raster.Rsave")

perennials_raster <- sum.twgd.trait(perennials,twgd_data,all_vars)
save(perennials_raster, file="global_map_life_form/perennials_raster.Rsave")

#########################
map_annuals <- GetSpRichness(annuals_raster)
save(map_annuals, file="global_map_life_form/map_annuals.Rsave")

map_perennials <- GetSpRichness(perennials_raster)
save(map_perennials, file="global_map_life_form/map_perennials.Rsave")


all_life_forms <- calc(stack(map_perennials, map_annuals), sum)
proportion_annuals <- map_annuals/all_life_forms
