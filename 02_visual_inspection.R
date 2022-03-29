# send names to gbif
setwd("~/Desktop/WCVP_special_issue/James_perennial_annual/life_history_houwie")
#rm(list=ls())
library(ape)
library(phytools)
library(data.table)

################
prepare.data<- function(trees, data, scoring, reference_table) {
  traits <- list()
  for(i in 1:length(trees)) {
    one_tree <- trees[[i]]
    one_label <- names(trees)[i]
    taxon_names <- one_tree$tip.label
    result_traits <- data.frame(species=taxon_names, life_form=NA)
    one_dataset <- subset(reference_table, reference_table$gbif_name %in% taxon_names)
    for(taxon_index in 1:length(taxon_names)) {
      if(any(reference_table$gbif_name==taxon_names[taxon_index])) {
        tmp_wcvp <- reference_table$wcvp_name[which(reference_table$gbif_name==taxon_names[taxon_index])]
        tmp_dataset <- subset(data, data$taxon_name%in%tmp_wcvp)
        life_form <- tail(names(sort(table(tmp_dataset$lifeform_description))),1)
        if(life_form=="") {
          # try second one
          life_form <- tail(names(sort(table(tmp_dataset$lifeform_description))),2)
          if(life_form[1]=="") {
            result_traits[taxon_index,2] <- "no_life_form_on_database"
          } else {
            result_traits[taxon_index,2] <- life_form[1] 
          }
        } else { 
          result_traits[taxon_index,2]  <- life_form
        }
      } else {
        result_traits[taxon_index,2] <- "no_life_form_on_database"
      }
      cat(taxon_index, "\r")
    }
    for(j in sequence(nrow(result_traits))) {
      life_form2 <- result_traits[j,2]
      if(life_form2!="no_life_form_on_database") {
        result_traits[j,2] <- scoring$scoring1[which(scoring$all_life_forms==life_form2)]        
      }
    }
    write.csv(result_traits, file=paste0("trait_dataset/",one_label,"_life_form.csv"), row.names = F)
    traits[[i]] <- result_traits
    names(traits)[i] <- names(trees)[i]
    
    plot.life.form(group_tree=one_tree, group_traits=result_traits, group=names(trees)[i])
    
    cat(one_label, " done.", "\n")
  }
  return(traits)
}

################
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
    pdf(paste0(getwd(), "/figures/", group,"_life_form_overview.pdf"), width= 4, height= 12)
    
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

################
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
################
#----------------
dist_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")

###################################
# Load trees
all_trees <- load.trees(tree.dir="trees")
focal_species_trees <- unname(unlist(lapply(all_trees, "[[", "tip.label")))
# Load taxize reference tables back
reference_table <- list.files("taxized_reference_tables", full.names = T)
reference_table <- do.call(rbind, lapply(reference_table, read.csv))
reference_table <- subset(reference_table, reference_table$gbif_name %in% focal_species_trees)
data <- subset(all_vars, all_vars$taxon_name %in% reference_table$wcvp_name)

###################################
scoring = read.csv("2022-03-19_life_form.csv")
trees = all_trees
# Scoring based on WCVP data
traits <- prepare.data(trees, data, scoring, reference_table) 

