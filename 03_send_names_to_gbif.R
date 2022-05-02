# send names to gbif
setwd("~/Desktop/WCVP_special_issue/James_perennial_annual/life_history_houwie")
library(ape)
library(rgbif)

get.tip.names <- function(tree_files) {
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
  return(focal_species_trees)
}

fix.names.taxize <- function(focal_species_trees) {
  for(name_index in 1:length(focal_species_trees)){
    one_tmp_string <- focal_species_trees[name_index]
    if(any(grepl("[()]", one_tmp_string))){
      splitted_names <- strsplit(one_tmp_string," ")[[1]]
      begin_author <- which(grepl("[()]", splitted_names))[1]
      species_name <- paste0(splitted_names[1:(begin_author-1)], collapse=" ")
      author <- splitted_names[begin_author:length(splitted_names)]
      old_authors <- author[grep("[()]", author)]
      end_first_half <- floor(length(old_authors)/2)
      before <- old_authors[1:end_first_half]
      after <- old_authors[(end_first_half+1):(length(old_authors))]
      if(paste(before,collapse = " ") == paste(after, collapse=" ")) {
        author <- paste(author[1:(length(author)/2)], collapse=" ")
        focal_species_trees[name_index] <- paste0(species_name, " ", author, collapse=" ")
      } else {
        author <- paste(author, collapse=" ")
        focal_species_trees[name_index] <- paste0(species_name, " ", author, collapse=" ")
      }
    }
  }
  return(focal_species_trees)
}

# Matching with trees
tree.dir <- "trees_gbif_tips"
tree_files <- list.files(tree.dir, full.names = T)
focal_species_trees <- get.tip.names(tree_files)
focal_species_trees <- fix.names.taxize(focal_species_trees)

test <- focal_species_trees[grep("Odontarrhena", focal_species_trees)]
test <- paste0(unlist(lapply(strsplit(test, " "), "[[", 1)), " ", unlist(lapply(strsplit(test, " "), "[[", 2)))
final <- resolveGBIF(test)

#focal_species_trees <- simplify.names.taxize(focal_species_trees)

# PILOT #----------------------
#focal_species_trees <- get.tip.names(tree_files[grep("Antirrhineae", tree_files)])
#----------------------

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


