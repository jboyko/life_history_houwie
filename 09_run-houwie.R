# imports
require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(data.table)

# functions
organizeData <- function(clade_name, climate_variable, data_files, tree_files){
  focal_tree_file <- tree_files[grep(clade_name, tree_files)]
  focal_data_file <- data_files[grep(clade_name, data_files)]
  focal_data_file <- focal_data_file[grep(paste0(climate_variable, "_"), focal_data_file)]
  phy <- read.tree(focal_tree_file)
  dat <- read.csv(focal_data_file)
  dat$species <- gsub(" ", "_", dat$species)
  dat <- dat[match(phy$tip.label, dat$species),]
  # cat("\n", "Removed", length(which(dat$life_form == "no_life_form_on_database")), "species of", length(dat$species), "because they didn't have life history data.\n")
  dat <- dat[!dat$life_form == "no_life_form_on_database",]
  plot_data <- data.frame(id = dat$species, 
                          life_form = as.factor(dat[,grep("life_form", colnames(dat))]),
                          value = dat[,grep("mean", colnames(dat))],
                          se = dat[,grep("se", colnames(dat))])
  # cat("\n", "Removed", length(which(apply(plot_data, 1, function(x) !any(is.na(x))))), "species of", length(dat$species), "because they didn't have climate data.\n")
  plot_data <- plot_data[apply(plot_data, 1, function(x) !any(is.na(x))),]
  pruned_phy <- keep.tip(phy, phy$tip.label[match(plot_data$id, phy$tip.label)])
  return(list(dat=plot_data, phy = pruned_phy))
}

runSingleModelSet <- function(clade_name, climate_variable, model_set, data_files, tree_files){
  focal_data <- organizeData(clade_name, climate_variable, data_files, tree_files)
  print(paste0("Running ", length(model_set), " models for ", clade_name, " (", length(focal_data$phy$tip.label), " taxa)."))
  model_set_res <- mclapply(model_set, function(x) hOUwie(focal_data$phy, focal_data$dat, ifelse(dim(x)[2] == 2, 1, 2), "ARD", x, nSim = 50, quiet = TRUE, mserr = "known"), mc.cores = 10)
  file_name <- paste0("res_files/", clade_name, "_", climate_variable, ".Rsave")
  save(model_set_res, file = file_name)
}

# run

# working directory
setwd("~/2022_life-history/")

data_files <- dir("datasets_final_for_hOUwie/full_datasets/", full.names = TRUE)
group_names <- unique(unlist(lapply(strsplit(dir("datasets_final_for_hOUwie/full_datasets/"), "-"), function(x) x[[1]])))

tree_files <- dir("trees_simplified_tips/", full.names = TRUE)
clade_name <- "Gesneriaceae"
# group_names[3]
# runSingleModelSet(clade_name, climate_variable, continuous_models[3], data_files, tree_files)

# climate variables
climatic_variables <- c(paste0("bio_", 1:19), "bio_ai", "bio_et0")

# continuous models
CID_model_names <- c("BM1", "OU1")
CD_model_names <- c("BMV", "OUA", "OUV", "OUM", "OUVA", "OUMV", "OUMA", "OUMVA")
CID2_model_names <- c("BMV", "OUA", "OUV", "OUM", "OUVA", "OUMV", "OUMA", "OUMVA")
CID_models <- lapply(CID_model_names, function(x) getOUParamStructure(x, 2, 1))
CD_models <- lapply(CD_model_names, function(x) getOUParamStructure(x, 2, 1))
CID2_models <- lapply(CID2_model_names, function(x) getOUParamStructure(x, 2, 2, TRUE))

names(CID_models) <- paste0("CID_", CID_model_names)
names(CD_models) <- paste0("CD_", CD_model_names)
names(CID2_models) <- paste0("CID_", CID2_model_names)

continuous_models <- c(CID_models, CD_models, CID2_models)

# run models
# bio5
climate_variable <- climatic_variables[5]
group_names_to_run <- group_names[!group_names %in% gsub("_.*", "", dir("res_files/"))]
mclapply(group_names_to_run, function(x) runSingleModelSet(x, climate_variable, continuous_models, data_files, tree_files), mc.cores = 5)


runSingleModelSet(group_names[3], climate_variable, continuous_models, data_files, tree_files)
hOUwie(focal_data$phy, focal_data$dat, 2, "ARD", continuous_models[[18]], nSim = 25, quiet = TRUE, diagn_msg = TRUE)




