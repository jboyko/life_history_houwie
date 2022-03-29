# 
setwd("~/Desktop/WCVP_special_issue/James_perennial_annual/life_history_houwie")
#rm(list=ls())

all_tables <- list.files("trait_datasets", full.names = T)
labels <- gsub(paste0(c("trait_datasets/","_life_form.csv"), collapse="|"), "", all_tables)

answer <- c()
for(i in 1:length(all_tables)) {
  one_dataset <- read.csv(all_tables[i])
  if(length(table(one_dataset$life_form))>2) {
    if(nrow(one_dataset)>50 && nrow(one_dataset)<2000) {
      answer <- c(answer, labels[i])
    }
  }
}

sink("good_trees_for_houwie.txt")
for(u in 1:length(answer)) {
  cat(answer[u],"\n")
}
sink()


# Extracting Balsamiaceae and Polemoniaceae from Ericales and Plantaginaceae from Lamiales
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

# crop Hyptidinae also for Andressa
all_trees <- load.trees(tree.dir="trees")
lamiales <- all_trees[[grep("Lamiales", names(all_trees))]]
ericales <- all_trees[[grep("Ericales", names(all_trees))]]

subset1 <- read.csv("cleaning_trees/polemoniaceae_tips.txt", h=F)
one_tree <- keep.tip(ericales, which(ericales$tip.label %in% as.character(subset1$V1)))
save(one_tree, file="Polemoniaceae-Rose_et_al-2018.Rsave")

subset2 <- read.csv("cleaning_trees/balsamiaceae_tips.txt", h=F)
one_tree <- keep.tip(ericales, which(ericales$tip.label %in% as.character(subset2$V1)))
save(one_tree, file="Balsamiaceae-Rose_et_al-2018.Rsave")


write.csv(lamiales$tip.label, "lamiales.csv")
write.csv(ericales$tip.label, "ericales.csv")



# (2) select trees that are okay for houwie 
# certain variation and number of tips

all_trees <- load.trees(tree.dir="trees")
good_trees <- read.table("good_trees_for_houwie.txt", h=F)
good_trees <- as.character(good_trees[,1])
trees <- subset(all_trees, names(all_trees)%in%good_trees)
reference_table <- list.files("taxized_reference_tables", full.names = T)
reference_table <- do.call(rbind, lapply(reference_table, read.csv))
reference_table <- subset(reference_table, reference_table$gbif_name %in% focal_species_trees)
data <- subset(all_vars, all_vars$taxon_name %in% reference_table$wcvp_name)
scoring = read.csv("2022-03-19_life_form.csv")

# Scoring based on WCVP data
traits <- prepare.data(trees, data, scoring, reference_table) 

