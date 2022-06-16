# rm(list=ls())
setwd("~/Desktop/WCVP_special_issue/James_perennial_annual/life_history_houwie")

library(data.table)
library(maptools)
library(raster)
library(sp)
library(rgeos)
library(rworldmap)
data("wrld_simpl")

###############################
Thinning <- function(points, species="scientificName", lat = "decimalLatitude", lon="decimalLongitude", n = 1) {
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  colnames(tmp_points)[colnames(tmp_points)==species] <- "tmp_sp"
  spp <- unique(tmp_points[,tmp_sp])
  results <- list()
  for(species_index in 1:length(spp)) {
    coords <- tmp_points[tmp_points[,tmp_sp]==spp[species_index],c("y","x")]
    coords <- coords[!duplicated(coords[,"x"]) & !duplicated(coords[,"y"]),]
    if(nrow(coords) > 1) {
      sp::coordinates(coords) <- ~ y + x
      raster::crs(coords) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
      r0 <- raster::raster(coords)
      raster::res(r0) <- 1 # cell resolution
      r0 <- raster::extend(r0, raster::extent(r0) + 5) 
      res <- cbind(spp[species_index], as.data.frame(dismo::gridSample(coords, r0, n))) # n = maximum number of points per cell
      colnames(res) <- c("tmp_sp", "lat","lon")
      results[[species_index]] <- res
    } else {
      res <- cbind(spp[species_index],coords)
      colnames(res) <- c("tmp_sp", "lat","lon")
      results[[species_index]] <- res
    }
  }
  results <- do.call(rbind, results)
  colnames(results) <- c(species, lat, lon)
  return(results)
}


###############################
GetClimateSummStats_custom <- function (points, type=c("raw","transformed")) {
  tmp_points <- points[,-which(colnames(points) %in% c("lon","lat"))]
  spp <- unique(tmp_points$species)
  vars <- colnames(tmp_points)[2]
  n_i<-c()
  sigma2_wi <- c()
  summ_stats <- matrix(nrow=length(spp), ncol=5)
    for(species_index in 1:length(spp)){
      sp1 <- tmp_points[tmp_points$species==spp[species_index],]
      cat("Now doing species", species_index, "\r")
      values <- sp1[,2]
      values <- values[!is.na(values)]
      if(type=="raw") {
        if(vars %in% c("bio_1","bio_4","bio_5","bio_6")) {
          values <-  (values / 10) 
        }
        n_i[species_index] <- length(values) # sample size
        sigma2_wi[species_index] <- ifelse(is.na(var(values)),0,var(values))  # sample variance
        
      }
      if(type=="transformed") {
        if(vars %in% c("bio_1","bio_4","bio_5","bio_6")) {
          values <-  (values / 10) + 273.15 # transforms to Kelvin
        }
        values <- log(values) # log
        n_i[species_index] <- length(values) # sample size
        sigma2_wi[species_index] <- ifelse(is.na(var(values)),0,var(values))  # sample variance
        
        if(any(values== -Inf)){
          values <- values[-which(values== -Inf)]
        }
      }
      n0 <- length(values)
      mean0 <- round(mean(values), 6)
      sd0 <- round(stats::sd(values), 6)
      se0 <- round(sd0/ sqrt(n0), 6)
      tmp_summ_stats <- c(n0, mean0, sd0, se0)
      summ_stats[species_index,] <- c(spp[species_index], tmp_summ_stats)
    }
  
  sigma2_w <- sum(sigma2_wi*(n_i - 1)) / sum(n_i - 1)
  within_sp_var <-  round(sigma2_w/n_i, 6)
  summ_stats <- cbind(summ_stats, within_sp_var)
  colnames(summ_stats)[6] <- paste0("within_sp_var_",vars)
  colnames(summ_stats) <- c("species",paste0("n_",vars), paste0("mean_",vars),
                            paste0("sd_",vars), paste0("se_",vars), paste0("within_sp_var_",vars))
  return(as.data.frame(summ_stats))
}

DataFromPoints <- function (points, layer) {
  if(any(colnames(points) != c("species","lat","lon"))) {
    stop("Columns have to be in the order of taxon, latitude and longitude and named as 'species', 'lat', and 'lon")
  }
  if(ncol(points)!=3) {
    stop("Dataset should be of class data.frame and organized in three columns named as 'species', 'lat', and 'lon'")   
  }
  if(!is.data.frame(points)) {
    stop("Dataset should be of class data.frame and organized in three columns named as 'species', 'lat', and 'lon'")   
  }
  cat("Extracting climatic information of", nrow(points), "points",  "\n")
  colnames(points) <- c("species", "lat", "lon")
  sp::coordinates(points) <- ~ lon + lat
  values <- raster::extract(layer, points)
  result <- cbind(points, values)
  out <- as.data.frame(result)
  if(class(layer@data@names) == "character"){
    colnames(out)[2] <- layer@data@names
  }
  return(out)
}

GetClimateSummStats <- function(climate_by_point, convert=TRUE){
  original_species_order <- unique(climate_by_point[,1])
  if(convert) {
    summ_stats <- aggregate(climate_by_point[,2], by = list(climate_by_point[,1]), function(x) BasicSummStats.k(x))
  } else {
    summ_stats <- aggregate(climate_by_point[,2], by = list(climate_by_point[,1]), function(x) BasicSummStats(x))
  }
  summ_stats <- cbind(summ_stats[,1], as.data.frame(summ_stats[,-1]))
  colnames(summ_stats) <- c("species", "n", "mean", "sd", "se")
  summ_stats <- summ_stats[match(original_species_order, summ_stats[,1]),]
  rownames(summ_stats) <- NULL
  return(summ_stats)
}

#BasicSummStats.old <- function(x){
#  x <- na.omit(x)
#  out <- c(n = length(x), 
#           mean = mean(x), 
#           sd = sd(x), 
#           se = sd(x)/sqrt(length(x)))
#  return(out)
#}


#-------------------------------
# Getting climate data
#-------------------------------
# Load cleaned points back:
all_cleaned_points <- fread("gbif_life_form/preliminary_cleaned_points.csv")

  # 1. Thinning occurence data first
thinned_points <- Thinning(all_cleaned_points, species="scientificName", lat = "decimalLatitude", lon="decimalLongitude", n = 1)

  # 2. Getting summary statistics of climatic variables for each species
colnames(thinned_points) <- c("species","lat","lon")

#########################################
all_layers <- list.files("climate_layers", ".tif$")
labels <- gsub(".tif$","", all_layers)
keep <- c("bio_1","bio_4","bio_5","bio_12","bio_13","bio_14","bio_15","bio_ai")
keep <- "bio_6"
all_layers <- subset(all_layers, labels %in% keep)
labels <- gsub(".tif$","", all_layers)
all_layers <- lapply(paste0("climate_layers/",all_layers), raster)
names(all_layers) <- labels

for(i in 1:length(all_layers)){
  one_layer <- all_layers[[i]]
  one_label <- names(all_layers)[i]
  allpoints <- DataFromPoints(thinned_points, one_layer)
  write.csv(allpoints, file=paste0("climate_data/",one_label,"_allpoints.csv"), row.names=F)
  summstats <- GetClimateSummStats_custom(allpoints, type="transformed")
  write.csv(summstats, file=paste0("climate_data/",one_label,"_summstats.csv"), row.names=F)
}



