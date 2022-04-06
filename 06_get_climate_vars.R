
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

#-------------------------------
# Getting climate data
#-------------------------------
# Load cleaned points back:
points.dir <- "./gbif_life_form"
cleaned_point_files <- list.files(points.dir, ".csv")
cleaned_point_files <- cleaned_point_files[grepl("cleaned", cleaned_point_files)]
all_cleaned_points <- lapply(paste0(points.dir, "/", cleaned_point_files), fread)
names(all_cleaned_points) <- unlist(lapply(strsplit(cleaned_point_files, "_"), "[[", 1))

# Directory to save preliminary datasets:
climate_data.dir <- "./final_datasets"
# Directory where climate layers are:
climate_layers.dir <- "./climate_layers"

{; for(family_index in 1:length(all_cleaned_points)) {
  # 1. Thinning occurence data first
  #thinned_points <- Thinning(all_cleaned_points[[family_index]], species="scientificName", lat = "decimalLatitude", lon="decimalLongitude", n = 1)
  # 2. Getting summary statistics of climatic variables for each species
  allpoints <- ClimateFromPoint_custom(thinned_points, species="scientificName",lon="decimalLongitude", lat="decimalLatitude", layerdir = climate_layers.dir)
  write.csv(allpoints, file=paste0(climate_data.dir, "/", names(all_cleaned_points)[family_index], "_allpoints.csv"), row.names=F)
  summstats <- GetClimateSummStats_custom(allpoints, type="raw")
  write.csv(summstats[[1]], file=paste0(climate_data.dir, "/", names(all_cleaned_points)[family_index], "_summstats_raw.csv"), row.names=F)
  summstats <- GetClimateSummStats_custom(allpoints, type="transformed")
  write.csv(summstats[[1]], file=paste0(climate_data.dir, "/", names(all_cleaned_points)[family_index], "_summstats.csv"), row.names=F)
}
  beepr::beep("fanfare"); } 

#-------------------------------
# Getting organized table for hOUwie
#-------------------------------
climate_data.dir <- "./final_datasets"
summstats_files <- list.files(climate_data.dir, "summstats.csv")
summstats <- lapply(paste0(climate_data.dir, "/", summstats_files), read.csv)
names(summstats) <- unlist(lapply(strsplit(summstats_files, "_"), "[[", 1))

# Trait data
trait.dir <- "./trait_dataset"
trait_files <- list.files(trait.dir, ".csv")
trait_files <- trait_files[grep("Antirrhineae", trait_files)]
traits <- lapply(paste0(trait.dir, "/", trait_files), read.csv)
labels <- unique(unlist(lapply(strsplit(trait_files, "-"), "[[", 1)))
names(traits) <- labels

{; for(family_index in 1:length(labels)) {
  group <- names(traits)[family_index]
  group_traits <- traits[[group]][,1:2]
  group_traits$species <- fix.names.taxize(group_traits$species)
  group_summstats <- summstats[[grep(group, names(summstats))]]
  
  # Matching datasets
  #group_summstats$species <- sub(" ","_", group_summstats$species)
  merged_table <- merge(group_summstats, group_traits, by="species", all=T)
  cleaned_table <- merged_table[,c("species","n_temp","mean_temp","se_temp","within_sp_var_temp","life_form")]
  write.csv(cleaned_table, file=paste0(climate_data.dir,"/",group, "_niche.csv"), row.names=F)
}
  beepr::beep("fanfare"); } 

