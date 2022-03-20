# send names to gbif
setwd("~/Desktop/James_perennial_annual/life_history_houwie")
library(ape)
library(rgbif)

# Matching with trees
tree.dir <- "trees"
tree_files <- list.files(tree.dir, full.names = T)
all_trees <- list()
for(i in 1:length(tree_files)) {
  load(tree_files[i])
  if(exists("one_tree")) {
    all_trees[[i]] <- one_tree
    names(all_trees)[i] <- gsub(paste0(c(paste0(tree.dir,"/"), ".Rsave"), collapse="|"),"", tree_files[i])
    rm("one_tree")
  }
}
focal_species_trees <- unname(unlist(lapply(all_trees, "[[", "tip.label")))

# Make sure the WCVP tables are in the same folder and load them
dist_sample <- read.csv("wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|")
names_sample <- read.csv("wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|")
# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")
# Load reference back
reference_table <- read.csv("reference_table.csv") # read the reference table that we created earlier again
reference_table <- subset(reference_table, focal_species_trees %in% reference_table$gbif_name)
# Write table for scoring
data <- subset(all_vars, all_vars$taxon_name %in% reference_table$wcvp_name)
all_life_forms <- data$lifeform_description
write.csv(as.data.frame(table(all_life_forms)), file="life_form.csv", row.names=F)






focal_species_wcvp <- all_vars %in% reference_table$wcvp_name 

# Now we send a request to GBIF to download the points for this list of species 

user <- "" # username
pwd <- "" # password
email <- "@gmail.com" # email


rgbif::occ_download(rgbif::pred_in("scientificName", taxized_names),
                    pred_in("basisOfRecord", 'PRESERVED_SPECIMEN'),
                    pred("hasCoordinate", TRUE),
                    format = "SIMPLE_CSV", user=user,pwd=pwd,email=email) # Sending request to GBIF

# After this step, log in your GBIF account and manually download the table with distribution points.
# (save the citation link too)


