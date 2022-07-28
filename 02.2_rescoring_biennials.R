# rm(list=ls())
setwd("~/Desktop/WCVP_special_issue/James_perennial_annual/life_history_houwie")

library(data.table)
library(maptools)
library(raster)
library(sp)
library(rgeos)
library(rworldmap)
data("wrld_simpl")

all_trait_files <- list.files("trait_dataset_post_curation/", full.names = T)
all_trait <- lapply(all_trait_files, read.csv)

scoring <- read.csv("2022-03-19_life_form.csv")

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
    write.csv(result_traits, file=paste0("trait_dataset_pre_curation/",one_label,"_life_form.csv"), row.names = F)
    traits[[i]] <- result_traits
    names(traits)[i] <- names(trees)[i]
    
    plot.life.form(group_tree=one_tree, group_traits=result_traits, group=names(trees)[i])
    
    cat(one_label, " done.", "\n")
  }
  return(traits)
}
