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
  dat <- dat[dat$species %in% phy$tip.label,]
  dat <- dat[match(phy$tip.label, dat$species),]
  # cat("\n", "Removed", length(which(dat$life_form == "no_life_form_on_database")), "species of", length(dat$species), "because they didn't have life history data.\n")
  dat <- dat[!dat$life_form == "no_life_form_on_database",]
  plot_data <- data.frame(id = dat$species, 
                          life_form = as.factor(dat[,grep("life_form", colnames(dat))]),
                          value = as.numeric(dat[,grep("mean", colnames(dat))]),
                          se = as.numeric(dat[,grep("se", colnames(dat))]))
  # cat("\n", "Removed", length(which(apply(plot_data, 1, function(x) !any(is.na(x))))), "species of", length(dat$species), "because they didn't have climate data.\n")
  plot_data <- plot_data[apply(plot_data, 1, function(x) !any(is.na(x))),]
  pruned_phy <- keep.tip(phy, phy$tip.label[phy$tip.label %in% plot_data$id])
  return(list(dat=plot_data, phy = pruned_phy))
}

# runSingleModelSet <- function(clade_name, climate_variable, model_set, data_files, tree_files){
#   print(paste0("Running ", length(model_set), " models for ", clade_name))
#   mclapply(names(model_set), function(x) runSingleModel(clade_name, climate_variable, x, model_set, data_files, tree_files), mc.cores = 10)
# }

runSingleModel <- function(clade_name, climate_variable, model_name, model_set, data_files, tree_files){
  focal_data <- organizeData(clade_name, climate_variable, data_files, tree_files)
  directory_path_a <- paste0("res_files/", climate_variable)
  dir.create(directory_path_a, showWarnings = FALSE)
  directory_path_b <- paste0("res_files/", climate_variable, "/", clade_name)
  dir.create(directory_path_b, showWarnings = FALSE)
  cont_model <- model_set[[match(model_name, names(model_set))]]
  res <- hOUwie(focal_data$phy, focal_data$dat, ifelse(dim(cont_model)[2] == 2, 1, 2), "ARD", cont_model, nSim = 50, quiet = TRUE, mserr = "known")
  file_name <- paste0(directory_path_b, "/", model_name, "_", clade_name, "_", climate_variable, ".Rsave")
  if(!is.null(res)){
    save(res, file = file_name)
  }
}

check_reruns <- function(climate_variable, clade_name, continuous_models){
  some_folders <- dir(paste0("res_files/", climate_variable), full.names = TRUE)
  model_names <- names(continuous_models)
  focal_folder <- dir(paste0("res_files/", climate_variable), full.names = TRUE)[grep(clade_name, dir(paste0("res_files/", climate_variable), full.names = TRUE))]
  failed_vector <- sapply(model_names, function(x) length(grep(paste0(x, "_"), dir(focal_folder))))
  if(any(failed_vector == 0)){
    df_rerun <- data.frame(clade = clade_name, model_rerun = names(failed_vector)[failed_vector == 0])
  }else{
    df_rerun <- c()
  }
  return(df_rerun)
}

complile_model_list <- function(climate_variable, focal_clade, model_names){
  print(focal_clade)
  some_folders <- dir(paste0("res_files/", climate_variable), full.names = TRUE)
  some_files <- dir(some_folders[grep(focal_clade, some_folders)], full.names = TRUE)
  complete_list <- list()
  for(i in 1:length(some_files)){
    load(some_files[i])
    if(is.null(res)){
      complete_list[[i]] <- NA
    }else{
      complete_list[[i]] <- res
    }
  }
  names(complete_list) <- unlist(lapply(strsplit(dir(some_folders[grep(focal_clade, some_folders)]), "_"), function(x) paste0(x[1], "_", x[2])))
  complete_list <- complete_list[match(model_names, names(complete_list))]
  complete_list <- complete_list[!unlist(lapply(complete_list, class)) == "logical"]
  directory_path <- paste0("res_files/compiled_models/", climate_variable)
  dir.create(directory_path, showWarnings = FALSE)
  file_name <- paste0(directory_path, "/", focal_clade, "_", climate_variable, ".Rsave")
  save(complete_list, file = file_name)
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
climate_variable <- "bio_5"
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

model_set <- continuous_models
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # run intial models # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
climate_variable <- "bio_1"
# clade_name <- "Alysseae"
models_to_run <- do.call(rbind, lapply(group_names, function(x) check_reruns(climate_variable, x, continuous_models)))
models_to_run <- do.call(c, apply(models_to_run, 1, list))
mclapply(models_to_run[1:10], function(x) runSingleModel(x[1], climate_variable, x[2], continuous_models, data_files, tree_files), mc.cores = 10)

# mclapply(group_names, function(x) runSingleModelSet(x, climate_variable, continuous_models, data_files, tree_files), mc.cores = 4)

# group_names_to_run <- group_names[!group_names %in% gsub("_.*", "", dir("res_files/"))]

# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # rerun  models didn't complete # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
climate_variable <- "bio_4"
models_to_rerun <- do.call(rbind, lapply(group_names, function(x) check_reruns(climate_variable, x, continuous_models)))
# models_to_rerun <- models_to_rerun[!models_to_rerun[,1] == group_names_to_run[1],]
# models_to_rerun <- models_to_rerun[!models_to_rerun[,1] == group_names_to_run[2],]
# models_to_rerun <- models_to_rerun[!models_to_rerun[,1] == group_names_to_run[3],]
# models_to_rerun <- models_to_rerun[!models_to_rerun[,1] == group_names_to_run[4],]
models_to_rerun <- do.call(c, apply(models_to_rerun, 1, list))

mclapply(models_to_rerun, function(x) runSingleModel(x[1], climate_variable, x[2], continuous_models, data_files, tree_files), mc.cores = 6)

# runSingleModel(models_to_rerun[[1]][1], climate_variable, models_to_rerun[[1]][2], continuous_models, data_files, tree_files)


# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # for running a missing clade # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
climate_variable <- "bio_5"
group_names_to_run <- group_names[!group_names %in% gsub("_.*", "", dir(paste0("res_files/",climate_variable)))]
models_to_rerun <- models_to_rerun[  models_to_rerun[,1] == group_names_to_run[1] |
                                     models_to_rerun[,1] == group_names_to_run[2] |
                                     models_to_rerun[,1] == group_names_to_run[3] |
                                     models_to_rerun[,1] == group_names_to_run[4],]
models_to_rerun <- do.call(c, apply(models_to_rerun, 1, list))

mclapply(models_to_rerun, function(x) runSingleModel(x[1], climate_variable, x[2], continuous_models, data_files, tree_files), mc.cores = 18)

# models_to_rerun <- models_to_rerun[models_to_rerun[,1] == group_names_to_run[2],]
# models_to_rerun <- models_to_rerun[models_to_rerun[,1] == group_names_to_run[3],]
# models_to_rerun <- models_to_rerun[models_to_rerun[,1] == group_names_to_run[4],]
# runSingleModelSet(group_names_to_run, climate_variable, continuous_models, data_files, tree_files)


# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # for unzipping the compiled models # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# climate_variable <- "bio_ai"
# some_files <- dir("res_files/compiled_models/", full.names = TRUE)
# some_files <- some_files[grep(paste0(climate_variable, ".Rsave"), some_files)]
# directory_path_a <- paste0("res_files/", climate_variable)
# dir.create(directory_path_a, showWarnings = FALSE)
# 
# for(i in 1:length(some_files)){
#   some_file <- some_files[i]
#   load(some_file)
#   print(some_file)
#   good_runs <- which(unlist(lapply(model_set_res, class)) == "houwie")
#   clade_name <- gsub("_.*", "", gsub(".*/", "", some_file))
#   directory_path_b <- paste0("res_files/", climate_variable, "/", clade_name)
#   dir.create(directory_path_b, showWarnings = FALSE)
#   for(j in good_runs){
#     model_name <- names(model_set_res)[j]
#     file_name <- paste0(directory_path_b, "/", model_name, "_", clade_name, "_", climate_variable, ".Rsave")
#     res <- model_set_res[[j]]
#     save(res, file = file_name)
#   }
# }

# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # for zipping the uncompiled models # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# bio_1  bio_14  bio_4  bio_ai bio_15 / bio_5 bio_12
climate_variable <- "bio_5"
# focal_clade <- "Lepidieae"
model_names <- names(continuous_models)
# complile_model_list(climate_variable, focal_clade, model_names)

lapply(group_names, function(x) complile_model_list(climate_variable, x, model_names))






