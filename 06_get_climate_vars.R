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

################################
ClimateFromPoint_custom <- function(points, species="scientificName",lon="lon", lat="lat", layerdir = ""){
  tmp_points = points
  colnames(tmp_points)[which(colnames(tmp_points) == lon)] <- "lon"
  colnames(tmp_points)[which(colnames(tmp_points) == lat)] <- "lat"
  colnames(tmp_points)[which(colnames(tmp_points) == species)] <- "species"
  tmp_points <- tmp_points[,c("species","lat","lon")]
  # Load climatic layers
  temp <- raster::raster(paste0(layerdir, "/bio_2.tif"))
  #prec <- raster::raster(paste0(layerdir, "/current_30sec/bio_12.tif"))
  #pet <- raster::raster(paste0(layerdir, "/et0_yr/et0_yr.tif"))
  #aridity <- raster::raster(paste0(layerdir, "/ai_et0/ai_et0.tif"))
  bio <- list(temp)
  vars <- c(temp)
  names(vars) <- c("temp")
  final_matrix <- matrix(nrow=nrow(tmp_points), ncol=length(vars))
  cat("Extracting climatic information of", nrow(tmp_points), "points",  "\n")
  sp::coordinates(tmp_points) <- ~ lon + lat
  for(var_index in 1:length(vars)) {
    layer <- bio[[var_index]]
    cat("\r",names(vars)[var_index])
    cat("","\n")
    values <- raster::extract(layer, tmp_points)
    final_matrix[,var_index] <- values
  }
  colnames(final_matrix) <- names(vars)
  result <- cbind(tmp_points, final_matrix)
  return(as.data.frame(result))
}

###############################
GetClimateSummStats_custom <- function (points, type=c("raw","transformed")) {
  tmp_points <- points[,-which(colnames(points) %in% c("lon","lat"))]
  vars <- c("temp")
  allclimatevars <- list()
  spp <- unique(tmp_points$species)
  for(var_index in 1:length(vars)) {
    cat("\r",vars[var_index])
    cat("","\n")
    n_i <- c()
    sigma2_wi <- c()
    summ_stats <- matrix(nrow=length(spp), ncol=5)
    for(species_index in 1:length(spp)){
      sp1 <- tmp_points[tmp_points$species==spp[species_index],]
      cat("\r","Now doing species", species_index)
      cat("","\n")
      values <- sp1[,vars[var_index]]
      values <- values[!is.na(values)]
      if(type=="raw") {
        if(vars[var_index] %in% c("temp")) {
          values <-  (values / 10) 
        }
        n_i[species_index] <- length(values) # sample size
        sigma2_wi[species_index] <- ifelse(is.na(var(values)),0,var(values))  # sample variance
        
      }
      if(type=="transformed") {
        if(vars[var_index] %in% c("temp")) {
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
      colnames(summ_stats) <- c("species",paste0("n_",vars[var_index]), paste0("mean_",vars[var_index]),
                                paste0("sd_",vars[var_index]), paste0("se_",vars[var_index]))
    }
    sigma2_w <- sum(sigma2_wi*(n_i - 1)) / sum(n_i - 1)
    within_sp_var <-  round(sigma2_w/n_i, 6)
    summ_stats <- cbind(summ_stats, within_sp_var)
    colnames(summ_stats)[6] <- paste0("within_sp_var_",vars[var_index])
    allclimatevars[[var_index]] <- summ_stats
  }
  return(allclimatevars)
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

GetClimateSummStats <- function(climate_by_point){
  original_species_order <- unique(climate_by_point[,1])
  summ_stats <- aggregate(climate_by_point[,2], by = list(climate_by_point[,1]), function(x) BasicSummStats(x))
  summ_stats <- cbind(summ_stats[,1], as.data.frame(summ_stats[,-1]))
  colnames(summ_stats) <- c("species", "n", "mean", "sd", "se")
  summ_stats <- summ_stats[match(original_species_order, summ_stats[,1]),]
  rownames(summ_stats) <- NULL
  return(summ_stats)
}

BasicSummStats <- function(x){
  x <- na.omit(x)
  out <- c(n = length(x), 
           mean = mean(x), 
           sd = sd(x), 
           se = sd(x)/sqrt(length(x)))
  return(out)
}
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
# Bio 15 ---------- add more layers later
layer <- raster("climate_layers/bio_15.tif")
allpoints <- DataFromPoints(thinned_points, layer)
write.csv(allpoints, file="climate_data/bio15_allpoints.csv", row.names=F)
summstats <- GetClimateSummStats(allpoints)
write.csv(summstats, file="climate_data/bio15_summstats.csv", row.names=F)
 
