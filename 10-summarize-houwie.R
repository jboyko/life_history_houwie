require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(data.table)
require(ggplot2)
require(reshape2)

# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # function # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
get_mod_avg_param_table <- function(climatic_variable, clade_name){
  print(clade_name)
  focal_files <- dir(res_folders[grep(paste0(climatic_variable, "$"), res_folders)], full.names = TRUE)
  focal_file <- focal_files[grep(clade_name, focal_files)]
  load(focal_file)
  mod_avg_param_table <- getModelAvgParams(complete_list, force = FALSE)
  return(mod_avg_param_table)
}

get_mod_table <- function(climatic_variable, clade_name){
  print(clade_name)
  focal_files <- dir(res_folders[grep(paste0(climatic_variable, "$"), res_folders)], full.names = TRUE)
  focal_file <- focal_files[grep(clade_name, focal_files)]
  load(focal_file)
  mod_avg_param_table <- getModelTable(complete_list)
  return(mod_avg_param_table)
}

getVariableName <- function(climatic_variable){
  bioclim <- c("BIO1 = Annual Mean Temperature",
               "BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))",
               "BIO3 = Isothermality (BIO2/BIO7) (×100)",
               "BIO4 = Temperature Seasonality (standard deviation ×100)",
               "BIO5 = Max Temperature of Warmest Month",
               "BIO6 = Min Temperature of Coldest Month",
               "BIO7 = Temperature Annual Range (BIO5-BIO6)",
               "BIO8 = Mean Temperature of Wettest Quarter",
               "BIO9 = Mean Temperature of Driest Quarter",
               "BIO10 = Mean Temperature of Warmest Quarter",
               "BIO11 = Mean Temperature of Coldest Quarter",
               "BIO12 = Annual Precipitation",
               "BIO13 = Precipitation of Wettest Month",
               "BIO14 = Precipitation of Driest Month",
               "BIO15 = Precipitation Seasonality (Coefficient of Variation)",
               "BIO16 = Precipitation of Wettest Quarter",
               "BIO17 = Precipitation of Driest Quarter",
               "BIO18 = Precipitation of Warmest Quarter",
               "BIO19 = Precipitation of Coldest Quarter",
               "BIOAI = Aridity Index")
  names(bioclim) <- c(paste0("bio_", 1:19), "bio_ai")
  return(bioclim[grep(climatic_variable, names(bioclim))])
}

# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # setup # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #

setwd("~/2022_life-history/")

res_folders <- dir("res_files/compiled_models/", full.names = TRUE)
group_names <- unique(unlist(lapply(strsplit(dir("datasets_final_for_hOUwie/full_datasets/"), "-"), function(x) x[[1]])))

data_files <- dir("datasets_final_for_hOUwie/full_datasets/", full.names = TRUE)

tree_files <- dir("trees_simplified_tips/", full.names = TRUE)

# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # run  # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #

climatic_variables <- c("bio_1", "bio_4", "bio_5", "bio_14", "bio_ai")

for(i in 1:length(climatic_variables)){
  climatic_variable <- climatic_variables[i]
  print("loading")
  mod_avg_param_list <- lapply(group_names, function(x) get_mod_avg_param_table(climatic_variable, x))
  names(mod_avg_param_list) <- group_names
  
  # get_mod_table("bio_4", "Balsamiaceae")
  
  # test <- mod_avg_param_list[[1]]
  # aggregate(test[,1:4], by = list(test$tip_state), mean)
  plot_data <- do.call(rbind, lapply(mod_avg_param_list, function(x) aggregate(x[,1:4], by = list(x$tip_state), mean)))
  plot_data$clade <- gsub("\\..*", "", rownames(plot_data))
  plot_data <- melt(plot_data, by = list("waiting_times", "alpha", "sigma.sq", "theta"))
  
  head(plot_data)
  
  title <- getVariableName(climatic_variable)
  print("plotting")
  ggplot(plot_data, aes(x = Group.1, y = value, group = clade, color = Group.1)) +
    ylab("") +
    xlab("Life history strategy") +
    ggtitle(title) +
    geom_line(color = "black") +
    geom_point(shape = 19) +
    facet_wrap(~variable, scales = "free") + 
    theme_bw()
  
  file_name <- paste0("figures/prelim_results/", climatic_variable, ".pdf")
  ggsave(file_name, height = 8, width = 14, units = "in")
}

# climatic_variable <- "bio_ai"
# clade_name <- "Solanaceae"
# focal_files <- dir(res_folders[grep(paste0(climatic_variable, "$"), res_folders)], full.names = TRUE)
# focal_file <- focal_files[grep(clade_name, focal_files)]
# load(focal_file)
# getModelTable(complete_list)



