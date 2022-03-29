# send names to gbif
setwd("~/Desktop/WCVP_special_issue/James_perennial_annual/life_history_houwie")
library(ape)
library(rgbif)

# Matching with trees
tree.dir <- "trees"
tree_files <- list.files(tree.dir, full.names = T)
#good_trees <- read.table("good_trees_for_houwie.txt", h=F)
#good_trees <- as.character(good_trees[,1])
#tree_files <- tree_files[gsub(paste0(c("trees/",".Rsave"), collapse="|"),"", tree_files) %in% good_trees]

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

# Now we send a request to GBIF to download the points for this list of species 
user <- "" # username
pwd <- "" # password
email <- "@gmail.com" # email

rgbif::occ_download(rgbif::pred_in("scientificName", focal_species_trees),
                    pred_in("basisOfRecord", 'PRESERVED_SPECIMEN'),
                    pred("hasCoordinate", TRUE),
                    format = "SIMPLE_CSV", user=user,pwd=pwd,email=email) # Sending request to GBIF

# After this step, log in your GBIF account and manually download the table with distribution points.
# (save the citation link too)


