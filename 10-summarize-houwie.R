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

get_stderr <- function(vec){
  return(sd(vec)/sqrt(length(vec)))
}

runTtest <- function(phy, plot_data, variable){
  annual_data <- plot_data[plot_data$Group.1 == "annual" & plot_data$variable == variable,]
  annual_vec <- annual_data$value
  names(annual_vec) <- annual_data$clade
  perrenial_data <- plot_data[plot_data$Group.1 == "perennial" & plot_data$variable == variable,]
  perrenial_vec <- perrenial_data$value
  names(perrenial_vec) <- perrenial_data$clade
  
  test <- phyl.pairedttest(phy, annual_vec, perrenial_vec, lambda = 1)
  return(test)
}

getSummTable <- function(big_list, climatic_variable){
  parameters <- c("waiting_times", "alpha", "sigma.sq", "theta")
  plot_data <- do.call(rbind, lapply(big_list[[match(climatic_variable, names(big_list))]], function(x) aggregate(x[,1:4], by = list(x$tip_state), mean)))
  plot_data$clade <- gsub("\\..*", "", rownames(plot_data))
  plot_data <- melt(plot_data, by = list("waiting_times", "alpha", "sigma.sq", "theta"))
  mean_data <- aggregate(plot_data[,4], by = list(plot_data$Group.1, plot_data$variable), mean)
  all_ttests <- lapply(parameters[c(1,1,3,4)], function(x) runTtest(phy, plot_data, x))
  stats <- do.call(rbind, lapply(all_ttests, function(x) c(phylo_mean=x$dbar, se=x$se, pvalue=x$P.dbar)))
  stats[2,] <- NA # alpha not different
  init_df <- data.frame(
    climatic_variable = climatic_variable,
    parameter = parameters,
    annual = mean_data$x[mean_data$Group.1=="annual"], 
    perennial = mean_data$x[mean_data$Group.1=="perennial"])
  out <- cbind(init_df, sig = ifelse(stats[,3] < 0.05, "*", "-"), stats)
  return(out)
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

climatic_variables <- c("bio_1", "bio_4", "bio_5", "bio_6", "bio_12", "bio_14", "bio_15", "bio_ai")
# climatic_variables <- c("bio_6")
big_list <- list()
print <- TRUE

for(i in 1:length(climatic_variables)){
  climatic_variable <- climatic_variables[i]
  print("loading")
  mod_avg_param_list <- lapply(group_names, function(x) get_mod_avg_param_table(climatic_variable, x))
  names(mod_avg_param_list) <- group_names
  big_list[[i]] <- mod_avg_param_list
  # get_mod_table("bio_4", "Balsamiaceae")
  
  # test <- mod_avg_param_list[[1]]
  # aggregate(test[,1:4], by = list(test$tip_state), mean)
  if(print){
    plot_data <- do.call(rbind, lapply(mod_avg_param_list, function(x) aggregate(x[,1:4], by = list(x$tip_state), mean)))
    data_se <- do.call(rbind, lapply(mod_avg_param_list, function(x) aggregate(x[,1:4], by = list(x$tip_state), get_stderr)))
    plot_data$clade <- gsub("\\..*", "", rownames(plot_data))
    plot_data <- melt(plot_data, by = list("waiting_times", "alpha", "sigma.sq", "theta"))
    
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
}

names(big_list) <- climatic_variables

for(i in 1:length(big_list)){
  print(i)
  focal_list <- big_list[[i]]
  tmp_table <- do.call(rbind, lapply(focal_list, function(x) aggregate(x[,1:4], by = list(x$tip_state), mean)))
  tmp_table <- cbind(clade=gsub("\\..*", "", rownames(tmp_table)), tmp_table)
  rownames(tmp_table) <- NULL
  file_name <- paste0("tables/climate_variable/", climatic_variables[i], ".csv")
  write.csv(tmp_table, file = file_name, row.names = FALSE)
  for(j in 1:length(focal_list)){
    focal_thing <- focal_list[[j]]
    file_name <- paste0("tables/tip_rates/", names(big_list)[i], "-", names(focal_list)[j], ".csv")
    write.csv(focal_thing, file = file_name, row.names = FALSE)
  }
}

# climatic_variable <- "bio_ai"
# clade_name <- "Solanaceae"
# focal_files <- dir(res_folders[grep(paste0(climatic_variable, "$"), res_folders)], full.names = TRUE)
# focal_file <- focal_files[grep(clade_name, focal_files)]
# load(focal_file)
# getModelTable(complete_list)

# annual_data <- plot_data[plot_data$Group.1 == "annual" & plot_data$variable == "theta",]
# annual_vec <- annual_data$value
# names(annual_vec) <- annual_data$clade
# perrenial_data <- plot_data[plot_data$Group.1 == "perennial" & plot_data$variable == "theta",]
# perrenial_vec <- perrenial_data$value
# names(perrenial_vec) <- perrenial_data$clade

phy <- read.tree("backbone_tree.tre")
phy$tip.label <- gsub("-.*", "", phy$tip.label)
climatic_variable <- "bio_1"


summ_tables <- do.call(rbind, lapply(climatic_variables, function(x) getSummTable(big_list, x)))
write.csv(summ_tables, file = "tables/summ_table.csv")

getSummTable(big_list, "bio_ai")









