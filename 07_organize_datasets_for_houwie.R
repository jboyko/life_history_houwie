#-------------------------------
# Getting organized table for hOUwie
#-------------------------------
# A way to simplify names in table and trees so that species names match again
simplify.names.taxize <- function(names) {
  results <- c()
  for(name_index in 1:length(names)){
    one_tmp_string <- names[name_index]
    splitted_names <- strsplit(one_tmp_string," ")[[1]]
    genus <- splitted_names[1]
    epiphet <- splitted_names[2]
    if(any(grepl("indet_sp",splitted_names))) {
      full_name <- "tip_to_drop" # indet species
    } else {
      if(stringr::str_detect(epiphet,"[[:upper:]]")) {
        full_name <- "tip_to_drop" # indet species
      } else {
        if(length(splitted_names) > 2) {
          complement <- splitted_names[3:length(splitted_names)]
          if(grepl("[()]", complement[1])) {
            full_name <- paste(c(genus, epiphet), collapse = " ")
          } else {
            if(stringr::str_detect(complement[1],"[[:upper:]]")) {
              full_name <- paste(c(genus, epiphet), collapse = " ")
            } else {
              complement <- subset(complement, !stringr::str_detect(complement,"[[:upper:]]"))
              complement <- subset(complement, !grepl(paste(c("[()]","&","([0-9]+).*$","^ex$"), collapse="|"), complement))
              if(length(complement)==0){
                full_name <- paste(c(genus, epiphet), collapse = " ")
              } else {
                full_name <- paste(c(genus, epiphet, complement), collapse = " ")
              }
            }
          } 
        }
      }
    }
    results[name_index] <- full_name
  }
  return(results)
}

# rm(list=ls())
setwd("~/Desktop/WCVP_special_issue/James_perennial_annual/life_history_houwie")
library(ape)

# Climate data
climate_data.dir <- "./climate_data"
summstats_files <- list.files(climate_data.dir, "summstats")
summstats <- lapply(paste0(climate_data.dir, "/", summstats_files), read.csv)
names(summstats) <- gsub("_summstats.csv","",summstats_files)
for(i in 1:length(summstats)) {
  one_var <- summstats[[i]] 
  one_var$species <- simplify.names.taxize(one_var$species)
  summstats[[i]] <- one_var
}

# Trait data
trait.dir <- "./trait_dataset_post_curation"
trait_files <- list.files(trait.dir, ".csv")
traits <- lapply(paste0(trait.dir, "/", trait_files), read.csv)
labels <- gsub("_life_form_cleaned.csv","",trait_files)
names(traits) <- labels

# Load simplified trees
tree.dir <- "./trees_simplified_tips/"
tree_files <- list.files(tree.dir, ".tre")
trees <- lapply(paste0(tree.dir, "/", tree_files), read.tree)
labels <- gsub("_cleaned.tre","",tree_files)
names(trees) <- labels

# Creating tables per groups
for(group_index in 1:length(trees)) {
  one_tree <- trees[[group_index]]
  one_label <- names(trees)[group_index]
  one_trait_dataset <- traits[[which(names(traits)==one_label)]]
  one_trait_dataset <- one_trait_dataset[,c("species","life_form")]
  one_trait_dataset$species <- gsub("_"," ",one_trait_dataset$species)
  for(var_index in 1:length(summstats)) {
    one_var <- summstats[[var_index]]
    merged <- merge(one_trait_dataset, one_var, by="species", all=T)
    merged <- subset(merged, !is.na(merged$life_form))
    write.csv(merged, file=paste0("datasets_final_for_hOUwie/full_datasets/", one_label,"_", names(summstats)[var_index],"_houwie_full.csv"), row.names = F)
  }
}

