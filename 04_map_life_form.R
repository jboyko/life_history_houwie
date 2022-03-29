# send names to gbif
setwd("~/Desktop/WCVP_special_issue/James_perennial_annual/life_history_houwie")
#rm(list=ls())
library(ape)
library(phytools)
library(data.table)
library(maptools)
library(raster)

dist_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")


reference_table <- list.files("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/taxized_reference_tables", full.names = T)
reference_table <- do.call(rbind, lapply(reference_table, read.csv))
scoring = read.csv("2022-03-19_life_form.csv")

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


path_tdwg = "/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools//wgsrpd-master/level3/level3.shp"
twgd_data <- suppressWarnings(maptools::readShapeSpatial(path_tdwg))

result_traits <- subset(result_traits, result_traits$life_form!="no_life_form_on_database")
  

occ_areas <- wcvp_subset$area_code_l3
    
    area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% occ_areas),]
    
    area_plus_buffer[1,][1]
    
  test <- raster(area_plus_buffer[1,])
  
  test[]  <- 1  
  over(area_plus_buffer[1,], test)
  plot(test)
  
  (area_plus_buffer[1,])
    point()
    
    class(area_plus_buffer)
    
    
    coords <- gbif_subset[,c("x","y")]
    sp::coordinates(coords) <- ~ x + y
    answer <- which(is.na(sp::over(coords, area_plus_buffer)[,3]))
    if(length(answer) != 0) {
      dubiousGBIF_ids <- c(dubiousGBIF_ids, as.character(gbif_subset$gbifID[answer]))
    }
  
  cat(species_index, "\r")

  
  
  library(sp)
  library(raster)
  
  ### Get Washington County map ###
  USCounty <- raster::getData('GADM', country='USA', level=2) # Download US County Map
  WA <- USCounty[USCounty$NAME_1 == "Washington",] # Subset data to just the State of WA
  
  ### Create point data ###
  
  # Set parameters
  numPts <- 30 # Set the number of points
  
  # Set projection I want
  ptsCRS <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
  
  # Project data
  WA <- spTransform(WA, ptsCRS)
  
  # Create point coordinates and attribute data
  set.seed(22)
  x <- sample(xmin(WA):xmax(WA), numPts, replace = TRUE) #xvalues in the extent of the state
  y <- sample(ymin(WA):ymax(WA), numPts, replace = TRUE) #yvalues in the extent of the state
  ptsCoords <- cbind(x,y) # Store the coordinates
  ptsData <- data.frame(Att = rnorm(numPts, 7, 1)) # Create and store attribute data
  
  pts <- SpatialPointsDataFrame(ptsCoords, data = ptsData, proj4string = crs(WA))
  
  # Plot the data
  plot(WA)
  plot(pts, add = TRUE)
  
  ### Perform Aggregation of point attributes by WA polygons ###
  countyAgg <- aggregate(pts[WA, "Att"], WA, sum)
  
  ### Plot Data ###
  spplot(countyAgg)

  assign


