# cleaning  points
# rm(list=ls())
setwd("~/Desktop/WCVP_special_issue/James_perennial_annual/life_history_houwie")

library(data.table)
library(maptools)
library(raster)
library(sp)
library(rgeos)
library(rworldmap)
data("wrld_simpl")

#-----------------------------
# If local
setwd("~/Desktop/WCVP_special_issue/James_perennial_annual/life_history_houwie/")
source("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/WCVPtools_functions.R")
source("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/fix_taxized_names.R")
dist_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
#-----------------------------
# If labcomputer
setwd("~/life_history_houwie")
source("../WCVPtools/WCVPtools_functions.R")
source("../WCVPtools/fix_taxized_names.R")
dist_sample <- read.table("../life_history_houwie/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("../life_history_houwie/wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
#-----------------------------
# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")

# reference table for taxized names
#-----------------------------
# If local
reference_table <- list.files("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/taxized_reference_tables", full.names = T)
reference_table <- do.call(rbind, lapply(reference_table, read.csv))
#-----------------------------
# If labcomputer
reference_table <- list.files("../WCVPtools/taxized_reference_tables", full.names = T)
reference_table <- do.call(rbind, lapply(reference_table, read.csv))

# Reading gbif file
gbif_data <- fread("gbif_life_form/all_points.csv") # load the table you downloaded from GBIF

gbif_data$scientificName[grep("Odontarrhena", gbif_data$scientificName)]
#gbif_data$scientificName[grep("Odontarrhena robertiana", gbif_data$scientificName)]

all_vars <- subset(all_vars, all_vars$genus %in% unique(gbif_data$genus))

# Looking at the WCVP table and TDWG to clean GBIF points
#-----------------------------
# If local
path="/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/wgsrpd-master/level3/level3.shp"
#-----------------------------
# If labcomputer
path="../WCVPtools/wgsrpd-master/level3/level3.shp"
#-----------------------------

twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))

cleaned_points <- gbif_data
#cleaned_points <- all_cleaned_points_files
#nrow(cleaned_points)
cleaned_points <- FilterWCVP_genus(cleaned_points, all_vars, twgd_data)
# testing a more "lose" cleaning
reference_table$gbif_name <- fix.names.taxize(reference_table$gbif_name)
subset_reference_table <- subset(reference_table, reference_table$gbif_name %in% unique(cleaned_points$scientificName))
if(nrow(subset_reference_table)>0){
  cleaned_points <- FilterWCVP(cleaned_points, all_vars, subset_reference_table, twgd_data) # This will filter the GBIF points acording to WCVP for species
}
# Cleaning common problems:
#cleaned_points <- RemoveNoDecimal(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
cleaned_points <- RemoveCentroids(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
cleaned_points <- RemoveDuplicates(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
cleaned_points <- RemoveZeros(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")

# Too much memomy required for this one, let's split the datasets:
cleaned_points1 <- RemoveSeaPoints(cleaned_points[1:floor(nrow(cleaned_points)/2),], lon="decimalLongitude", lat="decimalLatitude")
cleaned_points2 <- RemoveSeaPoints(cleaned_points[ceiling(nrow(cleaned_points)/2):nrow(cleaned_points),], lon="decimalLongitude", lat="decimalLatitude")
cleaned_points_final <- rbind(cleaned_points1, cleaned_points2)
write.csv(cleaned_points_final, file="gbif_life_form/preliminary_cleaned_points.csv", row.names=F)  

#------------------------
all_cleaned_points_files <- read.csv("gbif_life_form/preliminary_cleaned_points.csv")

# Plotting to inspect distributions
species <- unique(all_cleaned_points_files$genus)
species <- subset(species, species!="")
for(spp_index in 1:length(species)){
    tmp_subset <- as.data.frame(all_cleaned_points_files[all_cleaned_points_files$genus==species[spp_index],])
    pdf(paste0("distribution_maps/",species[spp_index],"_",unique(tmp_subset$family)[1],".pdf"))
    coord <- tmp_subset[,c("decimalLongitude","decimalLatitude")]
    coordinates(coord) <- ~ decimalLongitude + decimalLatitude
    plot(wrld_simpl)
    plot(coord, col="red", add=T)
    title(species[spp_index])
    cat(spp_index, "\r")
    dev.off()
  }


