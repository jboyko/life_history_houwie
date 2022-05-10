# Visual inspection of trees POST data curation (including both Kew's data and manually included data)
setwd("~/Desktop/WCVP_special_issue/James_perennial_annual/life_history_houwie")
#rm(list=ls())
library(ape)
library(phytools)
library(data.table)
library(phangorn)

load.trees <- function(tree.dir) {
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
  return(all_trees)
}

fix.names.taxize <- function(focal_species_trees) {
  for(name_index in 1:length(focal_species_trees)){
    one_tmp_string <- focal_species_trees[name_index]
    if(any(grepl("[()]", one_tmp_string))){
      splitted_names <- strsplit(one_tmp_string," ")[[1]]
      begin_author <- which(grepl("[()]", splitted_names))[1]
      species_name <- paste0(splitted_names[1:(begin_author-1)], collapse=" ")
      author <- splitted_names[begin_author:length(splitted_names)]
      author <- paste(author[1:(length(author)/2)], collapse=" ")
      focal_species_trees[name_index] <- paste0(species_name, " ", author, collapse=" ")
    } 
  }
  return(focal_species_trees)
}

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

# drop the negative branch lengths of the others
rm.br.length <- function(all_negbl, min_remain) {
  phy <- phy0 <- all_negbl
  while(min(phy$edge.length) < 0) {
    up_node <- phy$edge[which.min(phy$edge.length),][1]
    phy <- drop.tip(phy, Descendants(phy, up_node, type="tips")[[1]])
  }
  remain <- Ntip(phy) / Ntip(phy0)
  if(remain < min_remain) {
    return("phy should be dropped")
  } else {
   return(phy)
  }
}

plot.life.form <- function(group_tree, group_traits, group) {
  to_keep <- group_tree$tip.label[group_tree$tip.label %in% group_traits$species]
  if(length(to_keep)>1) {
    group_tree <- keep.tip(group_tree, to_keep)
    group_traits <- group_traits[group_traits$species %in% group_tree$tip.label,]
    group_traits <- group_traits[order(match(group_traits$species,group_tree$tip.label)),]
    
    mode <- group_traits$life_form
    names(mode) <- group_traits$species
    colors_states <- c()
    
    #tip.cols <- colors_states[as.factor(mode)]
    pdf(paste0(getwd(), "/figures/trees_post_curation/", group,"_life_form_overview.pdf"), width= 4, height= 12)
    
    plot(group_tree, show.tip.label=T, edge.width=0.2, adj=1, cex=0.08)
    par(fg="transparent")
    
    order <- colnames(to.matrix(mode, sort(unique(mode))))
    if(any("annual" %in% order)) {
      colors_states[which(order=="annual")] <- "midnightblue"
    }
    if(any("perennial" %in% order)) {
      colors_states[which(order=="perennial")] <- "goldenrod"
    }
    if(any("no_life_form_on_database" %in% order)) {
      colors_states[which(order=="no_life_form_on_database")] <- "grey"
    }
    tiplabels(pie=to.matrix(mode, sort(unique(mode))),piecol=colors_states,cex= 0.3,lwd=0.2, frame = "n")
    
    par(fg="black")
    #tiplabels(pch=1, bg=tip.cols, adj=1, cex=0.1, width = 0.1)
    
    legend("topleft", legend=sort(unique(mode)), pt.bg = colors_states, pch=21, cex=0.8)
    title(main=paste0(group))
    axisPhylo()
    dev.off()
  }
}

# Load trees
all_trees <- load.trees(tree.dir="trees_gbif_tips")

list_trees <- list()
# Create trees with simplified names
for(i in 1:length(all_trees)) {
  one_tree <- all_trees[[i]]
  one_label <- names(all_trees)[i]
  if(one_label=="Poaceae-Spriggs_et_al-2014"){
    next
  } else {
    one_tree$tip.label <- simplify.names.taxize(one_tree$tip.label)
    one_tree <- drop.tip(one_tree, which(one_tree$tip.label=="tip_to_drop"))
    if(min(one_tree$edge.length)<0){
      one_tree <- rm.br.length(one_tree, min_remain = 0.8)
    }
    list_trees[[i]] <- one_tree
    names(list_trees)[i] <- one_label
    write.tree(one_tree, file=paste0("trees_simplified_tips/", one_label, "_cleaned.tre")) # saving trees post-curation    
  }
}


# Make life form datasets that match the trees
trait_data_files <- list.files("trait_dataset_pre_curation")
labels <- gsub(paste0(c("_life_form_curated.csv","_life_form.csv"), collapse = "|"), "", trait_data_files)
trait_data <- lapply(paste0("trait_dataset_pre_curation/", trait_data_files), read.csv)
names(trait_data) <- labels

list_traits <- list()
for(i in 1:length(trait_data)) {
  one_dataset <- trait_data[[i]]
  one_label <- names(trait_data)[i]
  one_dataset$species <- simplify.names.taxize(one_dataset$species)
  one_dataset <- subset(one_dataset, one_dataset$species!="tip_to_drop")
  list_traits[[i]] <- one_dataset
  names(list_traits)[i] <- one_label
  write.csv(one_dataset, file=paste0("trait_dataset_post_curation/", one_label, "_life_form_cleaned.csv"), row.names = F)
}

for(i in 1:length(list_trees)) {
  one_tree <- list_trees[[i]]
  one_label <- names(list_trees)[i]
  one_trait_dataset <- list_traits[[grep(one_label, names(list_traits))]]
  plot.life.form(one_tree, one_trait_dataset, one_label)
}

# Table with proportion of NAs etc
# Provide a table with: phylogeny, number of species sampled, sf, proportion of annuals, proportion of perennials, proportion of missing data

#treebank <- read.csv("cleaning_trees/treebank_info_simplified.csv")
#treebank <- subset(treebank, treebank$label%in%names(list_traits))
#treebank <- treebank[,c("label","sf_ingroup")]
#write.csv(treebank, file="trees_sfs.csv", row.names = F)

# Build a table with data description:
trees_sfs <- read.csv("data_description/trees_sfs.csv")
results_table <- as.data.frame(matrix(nrow=length(list_traits), ncol=6))
colnames(results_table) <- c("tree","ntip","sf","annuals","perennials","missing_data")
for(i in 1:length(list_traits)) {
  one_dataset <- list_traits[[i]]
  results_table$tree[i] <- names(list_traits)[i]
  results_table$ntip[i] <- length(one_dataset$life_form)
  one_sf <- trees_sfs$sf_ingroup[trees_sfs$label==results_table$tree[i]]
  if(length(one_sf)!=0){
      results_table$sf[i] <- one_sf
  }
  results_table$annuals[i] <- round(length(which(one_dataset$life_form=="annual")) / (results_table$ntip[i]), 2)
  results_table$perennials[i] <- round(length(which(one_dataset$life_form=="perennial")) / (results_table$ntip[i]), 2)
  results_table$missing_data[i] <- round(length(which(one_dataset$life_form=="no_life_form_on_database")) / (results_table$ntip[i]), 2)
}

write.csv(results_table, file="data_description/summary_table_life_form.csv", row.names = F)


