# send names to gbif
setwd("~/Desktop/James_perennial_annual/life_history_houwie")
#rm(list=ls())
library(ape)
library(phytools)
library(data.table)

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
dist_sample <- fread("wcvp_names_and_distribution_special_edition_2022 (1)/wcvp_distribution.txt", sep="|")
names_sample <- fread("wcvp_names_and_distribution_special_edition_2022 (1)/wcvp_names.txt")

sort(unique(names_sample$family))

# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")
# Load reference back
reference_table <- read.csv("reference_table.csv") # read the reference table that we created earlier again
reference_table <- subset(reference_table, focal_species_trees %in% reference_table$gbif_name)
# Write table for scoring
data <- subset(all_vars, all_vars$taxon_name %in% reference_table$wcvp_name)
# all_life_forms <- data$lifeform_description
# write.csv(as.data.frame(table(all_life_forms)), file=paste0(Sys.Date(),"_life_form.csv"), row.names=F)


###################################
# 
###################################
scoring <- read.csv("life_form.csv")
trees <- all_trees


prepare.data<- function(trees, data, scoring) {
  traits <- list()
  for(i in 1:length(trees)) {
    cat(i, "\r")
    one_tree<- trees[[i]]
    taxon_names <- one_tree$tip.label
    result_traits <- data.frame(species=taxon_names, life_form=NA)
    one_dataset <- subset(reference_table, reference_table$gbif_name %in% taxon_names)
    if(nrow(one_dataset)==0){
      result_traits[,2] <- "no_life_form_on_database"
    }
    if(nrow(one_dataset)>0){
      # make trait dataset usable
      for(taxon_index in 1:length(taxon_names)){
        tmp_dataset <- subset(data, data$taxon_name==taxon_names[taxon_index])
        
        result_traits[taxon_index,2] 
        
        for(i in 1:nrow(data)) {
          data$gbif_ref[i] <- reference_table$gbif_name[which(data$taxon_name[i] == reference_table$wcvp_name)]
          cat(i, "\r")
        }
        
        life_form <- tail(names(sort(table(tmp_dataset$lifeform_description))),1)
        if(life_form=="") {
          # try second one
          life_form <- tail(names(sort(table(tmp_dataset$lifeform_description))),2)
          if(life_form=="") {
            life_form <- "no_life_form_on_database"
            result_traits[taxon_index,1] <- taxon_names[taxon_index]
            result_traits[taxon_index,2] <- life_form
          }
        } else {
          result_traits[taxon_index,1] <- taxon_names[taxon_index]
          result_traits[taxon_index,2] <- scoring$scoring1[which(scoring$all_life_forms==life_form)]        
        }
      }
      colnames(result_traits) <- c("species","life_form")
      traits[[i]] <- result_traits
      names(traits)[i] <- names(trees)[i]
    }
  }
  trees[unlist(lapply(traits, is.null))] <- NULL
  traits[unlist(lapply(traits, is.null))] <- NULL
}


for(group_index in 1:length(traits)) {
  one_dataset <- traits[[group_index]]
  for(i in 1:nrow(one_dataset)) {
    one_dataset$gbif_name[i] <- reference_table$gbif_name[reference_table$wcvp_name 
                                                          == one_dataset$species[i]]
  }
  traits[[group_index]] <- one_dataset
}


# Tree plots
pdf(paste0(getwd(), "/figures/", group, "_life_form_overview.pdf"), width= 4, height= 12)
for(group_index in 1:length(trees)) {
  group <- names(traits)[group_index]
  group_traits <- traits[[group]]
  group_tree <- trees[[grep(group, names(trees))]]
  to_keep <- group_tree$tip.label[group_tree$tip.label %in% group_traits$gbif_name]
  if(length(to_keep)>1) {
  group_tree <- keep.tip(group_tree, to_keep)
  group_traits <- group_traits[group_traits$gbif_name %in% group_tree$tip.label,]
  group_traits <- group_traits[order(match(group_traits$gbif_name,group_tree$tip.label)),]
  
  mode <- group_traits$life_form
  names(mode) <- group_traits$gbif_name
  colors_states <- c("midnightblue", "goldenrod", "red", "green")
  
  #tip.cols <- colors_states[as.factor(mode)]
  
  plot(group_tree, show.tip.label=T, edge.width=0.2, adj=1, cex=0.08)
  par(fg="transparent")
  tiplabels(pie=to.matrix(mode, sort(unique(mode))),piecol=colors_states,cex= 0.3,lwd=0.2, frame = "n")
  par(fg="black")
  #tiplabels(pch=1, bg=tip.cols, adj=1, cex=0.1, width = 0.1)
  
  legend("topleft", legend=sort(unique(mode)), pt.bg = colors_states, pch=21, cex=0.8)
  title(main=paste0(group))
  axisPhylo()
  }
}
dev.off()


