# which groups to keep:


# functions
organizeData <- function(clade_name, data_files, tree_files){
  focal_tree_file <- tree_files[grep(clade_name, tree_files)]
  focal_data_file <- data_files[grep(clade_name, data_files)]
  phy <- read.tree(focal_tree_file)
  dat <- read.csv(focal_data_file)
  dat$species <- gsub(" ", "_", dat$species)
  dat <- dat[match(phy$tip.label, dat$species),]
  #cat("\n", "Removed", length(which(dat$life_form == "no_life_form_on_database")), "species of", length(dat$species), "because they didn't have life history data.\n")
  dat <- dat[!dat$life_form == "no_life_form_on_database",]
  plot_data <- data.frame(id = dat$species, value = dat[,"mean"], life_form = as.factor(dat[,"life_form"]))
  #cat("\n", "Removed", length(which(apply(plot_data, 1, function(x) !any(is.na(x))))), "species of", length(dat$species), "because they didn't have climate data.\n")
  plot_data <- plot_data[apply(plot_data, 1, function(x) !any(is.na(x))),]
  pruned_phy <- keep.tip(phy, phy$tip.label[match(plot_data$id, phy$tip.label)])
  prop_kept <- round(Ntip(pruned_phy) / Ntip(phy), 2)
  if( prop_kept < 0.5 ) {
    cat("\n", clade_name, "should perhaps be removed. Only",prop_kept,  "of the original dataset was kept after filterying.\n")
    return(clade_name)
  } else {
    return("ok")
  }
  #return(list(dat=plot_data[,c(1,3,2)], phy = pruned_phy))
}


data_files <- dir("datasets_final_for_hOUwie/full_datasets", full.names = TRUE)
group_names <- unlist(lapply(strsplit(dir("datasets_final_for_hOUwie/full_datasets/"), "-"), function(x) x[[1]]))
tree_files <- dir("trees_simplified_tips/", full.names = TRUE)

answer <- c()
for(i in 1:length(unique(group_names))) {
  clade_name <- unique(group_names)[i]
  tmp <- data_files[grep(clade_name, data_files)][1]
  answer[i] <- organizeData(clade_name, tmp, tree_files)  
}
answer <- subset(answer, answer != "ok")



