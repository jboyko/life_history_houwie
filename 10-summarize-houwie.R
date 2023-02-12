# rm(list=ls())

require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(data.table)
require(ggplot2)
require(reshape2)
require(ggplotify)
require(gridExtra)
require(ggtree)
require(aplot)
require(dpl)

# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # function # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
get_mod_avg_param_table <- function(climatic_variable, clade_name, type){
  print(clade_name)
  focal_files <- dir(res_folders[grep(paste0(climatic_variable, "$"), res_folders)], full.names = TRUE)
  focal_file <- focal_files[grep(clade_name, focal_files)]
  load(focal_file)
  mod_table <- get_mod_table(climatic_variable, clade_name, type=type)
  new_list <- complete_list
  # new_list <- complete_list[mod_table[,6] < 4]
  # if(length(new_list) == 1){
  #   dummy_model <- complete_list[[which.max(mod_table[,6])]]
  #   dummy_model$AIC <- dummy_model$AICc <- dummy_model$BIC <- Inf
  #   new_list <- c(new_list, list(dummy_model))
  # }
  mod_avg_param_table <- getModelAvgParams(new_list, type=type, force = TRUE)
  return(mod_avg_param_table)
}

load_file <- function(climatic_variable, clade_name){
  print(clade_name)
  focal_files <- dir(res_folders[grep(paste0(climatic_variable, "$"), res_folders)], full.names = TRUE)
  focal_file <- focal_files[grep(clade_name, focal_files)]
  load(focal_file)
  return(complete_list)
}


get_mod_table <- function(climatic_variable, clade_name, type=type){
  print(clade_name)
  focal_files <- dir(res_folders[grep(paste0(climatic_variable, "$"), res_folders)], full.names = TRUE)
  focal_file <- focal_files[grep(clade_name, focal_files)]
  load(focal_file)
  mod_avg_param_table <- getModelTable(complete_list, type=type)
  return(mod_avg_param_table)
}

convertVariable <- function(climatic_variable, exp_vector){
  if(climatic_variable == 'bio_1'){
    return(exp(exp_vector) - 273.15)
  }
  if(climatic_variable == 'bio_4'){
    return(exp(exp_vector))
  }
  if(climatic_variable == 'bio_5'){
    return(exp(exp_vector) - 273.15)
  }
  if(climatic_variable == 'bio_6'){
    return(exp(exp_vector) - 273.15)
  }
  if(climatic_variable == 'bio_12'){
    return(exp(exp_vector))
  }
  if(climatic_variable == 'bio_14'){
    return(exp(exp_vector))
  }
  if(climatic_variable == 'bio_15'){
    return(exp(exp_vector))
  }
  if(climatic_variable == 'bio_ai'){
    return(exp(exp_vector) * 0.0001)
  }
}

getVariableName <- function(climatic_variable){
  # bioclim <- c("BIO1 = Annual Mean Temperature",
  #              "BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))",
  #              "BIO3 = Isothermality (BIO2/BIO7) (×100)",
  #              "BIO4 = Temperature Seasonality (standard deviation ×100)",
  #              "BIO5 = Max Temperature of Warmest Month",
  #              "BIO6 = Min Temperature of Coldest Month",
  #              "BIO7 = Temperature Annual Range (BIO5-BIO6)",
  #              "BIO8 = Mean Temperature of Wettest Quarter",
  #              "BIO9 = Mean Temperature of Driest Quarter",
  #              "BIO10 = Mean Temperature of Warmest Quarter",
  #              "BIO11 = Mean Temperature of Coldest Quarter",
  #              "BIO12 = Annual Precipitation",
  #              "BIO13 = Precipitation of Wettest Month",
  #              "BIO14 = Precipitation of Driest Month",
  #              "BIO15 = Precipitation Seasonality (Coefficient of Variation)",
  #              "BIO16 = Precipitation of Wettest Quarter",
  #              "BIO17 = Precipitation of Driest Quarter",
  #              "BIO18 = Precipitation of Warmest Quarter",
  #              "BIO19 = Precipitation of Coldest Quarter",
  #              "BIOAI = Aridity Index")
  bioclim <- c("Annual Mean Temperature (°C)",
               "Mean Diurnal Range",
               "Isothermality (BIO2/BIO7) (×100)",
               "Temperature Seasonality (SD x 100)",
               "Max Temperature of Warmest Month (°C)",
               "Min Temperature of Coldest Month (°C)",
               "Temperature Annual Range (BIO5-BIO6)",
               "Mean Temperature of Wettest Quarter",
               "Mean Temperature of Driest Quarter",
               "Mean Temperature of Warmest Quarter",
               "Mean Temperature of Coldest Quarter",
               "Annual Precipitation (mm)",
               "Precipitation of Wettest Month (mm)",
               "Precipitation of Driest Month (mm)",
               "Precipitation Seasonality (Coefficient of Variation)",
               "Precipitation of Wettest Quarter",
               "Precipitation of Driest Quarter",
               "Precipitation of Warmest Quarter",
               "Precipitation of Coldest Quarter",
               "Aridity Index (P/PET)")
  names(bioclim) <- c(paste0("bio_", 1:19), "bio_ai")
  
  return(bioclim[match(climatic_variable, names(bioclim))])
}

get_stderr <- function(vec){
  return(sd(vec)/sqrt(length(vec)))
}

runTtest <- function(phy, plot_data, variable=3){
  annual_data <- plot_data[plot_data$Group.2 == "annual",]
  annual_vec <- annual_data[,variable]
  names(annual_vec) <- annual_data$Group.1
  perrenial_data <- plot_data[plot_data$Group.2 == "perennial",]
  perrenial_vec <- perrenial_data[,variable]
  names(perrenial_vec) <- perrenial_data$Group.1
  
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

makePlot <- function(variable, letter, mu, TABLE = FALSE){
  climate_variable <- getVariableName(variable)
  
  focal_list <- all_model_tables[[variable]]
  focal_dat <- do.call(rbind, focal_list)
  taxon_names <- unlist(lapply(focal_list, function(x) rownames(x)))
  clade_names <- gsub("\\..*", "", rownames(focal_dat))
  rownames(focal_dat) <- NULL
  focal_dat$taxon <- taxon_names
  focal_dat$clade <- clade_names
  # head(focal_dat)
  summ_data <- aggregate(focal_dat[,c(1:6)], by = list(focal_dat$clade, focal_dat$tip_state), mean)
  summ_data[,7] <- convertVariable(variable, summ_data[,7])
  summ_data[,6] <- convertVariable(variable, summ_data[,6])
  summ_data <- summ_data[!summ_data$Group.1 == "Chorisporeae",] # remove Chorisporeae
  write.csv(summ_data, paste0("tables/parameter_tables/", variable, ".csv"))
  if(TABLE){
    all_ttests <- lapply(c(3,5:8), function(x) runTtest(phy, summ_data, x))
    stats <- do.call(rbind, lapply(all_ttests, function(x) c(phylo_mean=x$dbar, se=x$se, pvalue=x$P.dbar)))
    out <- data.frame(climate_variable = as.character(climate_variable), model_variable = c("rate", "sigma.sq", "theta", "expected_mean", "expected_variance"), stats)
    return(out)
  }
  if(mu == "mean" | mu == "both"){
    ttest_res <- runTtest(phy, summ_data, 7)
    a <- ggplot(summ_data %>% group_by(Group.1) %>% mutate(slope = (expected_mean[Group.2=="annual"] - expected_mean[Group.2=="perennial"])/(2-1)), 
    aes(x = Group.2, y = expected_mean, group = Group.1, color = slope > 0)) +
      xlab("") +
      ylab("Expected mean") +
      geom_line() +
      geom_point(shape = 19) +
      ggtitle(paste0(letter, ") ", climate_variable), subtitle = paste0("p = ", round(ttest_res$P.dbar, 3))) +
      theme_bw() +
      stat_summary(fun=mean, geom="point",aes(group=1, size = 2)) +  
      stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group=1), width = 0.15, color = "black") +
      # labs(caption = paste0("p = ", round(t_sigma$P.dbar, 3))) +
      theme(legend.position="none", text = element_text(size = 15), axis.text.y = element_text(size = 10)) +
      theme(plot.title=element_text(size=10)) +
      theme(plot.subtitle=element_text(size=10)) +
      theme(axis.title=element_text(size=10)) +
      theme(axis.text=element_text(size=10)) 
    out <- a
  }
  if(mu == "var" | mu == "both"){
    ttest_res <- runTtest(phy, summ_data, 8)
    b <- ggplot(summ_data %>% group_by(Group.1) %>% mutate(slope = (expected_var[Group.2=="annual"] - expected_var[Group.2=="perennial"])/(2-1)), 
    aes(x = Group.2, y = expected_var, group = Group.1, color = slope > 0)) +
      xlab("") +
      ylab("Expected variance") +
      geom_line() +
      geom_point(shape = 19) +
      ggtitle(paste0(letter, ") ", climate_variable), subtitle = paste0("p = ", round(ttest_res$P.dbar, 3))) +
      theme(plot.title = element_text(size = 0.2)) +
      theme_bw() +
      stat_summary(fun=mean, geom="point",aes(group=1, size = 2)) +  
      stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group=1), width = 0.15, color = "black") +
      # labs(caption = paste0("p = ", round(t_sigma$P.dbar, 3))) +
      theme(legend.position="none", text = element_text(size = 15), axis.text.y = element_text(size = 10)) +
      theme(plot.title=element_text(size=10)) +
      theme(plot.subtitle=element_text(size=10)) +
      theme(axis.title=element_text(size=10)) +
      theme(axis.text=element_text(size=10)) 
    out <- b
  }
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
# # # # # initial summarization # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #

climatic_variables <- c("bio_1", "bio_4", "bio_5", "bio_6", "bio_12", "bio_14", "bio_15", "bio_ai")
# climatic_variables <- c("bio_6")
big_list <- list()
print <- FALSE

for(i in 1:length(climatic_variables)){
  climatic_variable <- climatic_variables[i]
  print("loading")
  # mod_avg_param_list <- lapply(group_names, function(x) get_mod_avg_param_table(climatic_variable, x, type="AICc"))
  # test <- load_file(climatic_variable, "Balsamiaceae")
  # mod_avg_param_list <- lapply(group_names, function(x) get_mod_table(climatic_variable, x))
  mod_avg_table_list <- lapply(group_names, function(x) get_mod_table(climatic_variable, x, "AICc"))
  names(mod_avg_table_list) <- group_names
  big_list[[i]] <- mod_avg_table_list
  # get_mod_table(climatic_variable, "CES")
  # test <- mod_avg_param_list[[1]]
  # mod_avg_param_list <- big_list[[i]]
  # aggregate(test[,1:4], by = list(test$tip_state), mean)
  # if(print){
  #   plot_data <- do.call(rbind, lapply(mod_avg_param_list, function(x) aggregate(x[,1:4], by = list(x$tip_state), mean)))
  #   data_se <- do.call(rbind, lapply(mod_avg_param_list, function(x) aggregate(x[,1:4], by = list(x$tip_state), get_stderr)))
  #   plot_data$clade <- gsub("\\..*", "", rownames(plot_data))
  #   title <- getVariableName(climatic_variable)
  #
  #   # discrete rate
  #   plot_data_1 <- plot_data[,c(1,2,6)]
  #   ranges <- c(mean(plot_data_1[,2]) - (sd(plot_data_1[,2])), mean(plot_data_1[,2]) + (sd(plot_data_1[,2])))
  #   p1 <- ggplot(melt(plot_data_1), aes(x = Group.1, y = value, group = clade, color = Group.1)) +
  #     ylab("") +
  #     xlab("Life history strategy") +
  #     ggtitle("Discrete rate") +
  #     geom_line(color = "black") +
  #     geom_point(shape = 19) +
  #     coord_cartesian(ylim=ranges) +
  #     theme_bw()
  #
  #   # alpha
  #   plot_data_1 <- plot_data[,c(1,3,6)]
  #   ranges <- c(mean(plot_data_1[,2]) - (sd(plot_data_1[,2])), mean(plot_data_1[,2]) + (sd(plot_data_1[,2])))
  #   p2 <- ggplot(melt(plot_data_1), aes(x = Group.1, y = value, group = clade, color = Group.1)) +
  #     ylab("") +
  #     xlab("Life history strategy") +
  #     ggtitle("Alpha") +
  #     geom_line(color = "black") +
  #     geom_point(shape = 19) +
  #     coord_cartesian(ylim=ranges) +
  #     theme_bw()
  #
  #   # sigma.sq
  #   plot_data_1 <- plot_data[,c(1,4,6)]
  #   ranges <- c(mean(plot_data_1[,2]) - (sd(plot_data_1[,2])), mean(plot_data_1[,2]) + (sd(plot_data_1[,2])))
  #   p3 <- ggplot(melt(plot_data_1), aes(x = Group.1, y = value, group = clade, color = Group.1)) +
  #     ylab("") +
  #     xlab("Life history strategy") +
  #     ggtitle("Sigma Squared") +
  #     geom_line(color = "black") +
  #     geom_point(shape = 19) +
  #     coord_cartesian(ylim=ranges) +
  #     theme_bw()
  #
  #   # sigma.sq
  #   plot_data_1 <- plot_data[,c(1,5,6)]
  #   ranges <- c(mean(plot_data_1[,2]) - (sd(plot_data_1[,2])), mean(plot_data_1[,2]) + (sd(plot_data_1[,2])))
  #   p4 <- ggplot(melt(plot_data_1), aes(x = Group.1, y = value, group = clade, color = Group.1)) +
  #     ylab("") +
  #     xlab("Life history strategy") +
  #     ggtitle("Theta") +
  #     geom_line(color = "black") +
  #     geom_point(shape = 19) +
  #     coord_cartesian(ylim=ranges) +
  #     theme_bw()
  #
  #   # plot_data <- melt(plot_data, by = list("waiting_times", "alpha", "sigma.sq", "theta"))
  #
  #   title <- getVariableName(climatic_variable)
  #   print("plotting")
  #   ggplot(plot_data, aes(x = Group.1, y = value, group = clade, color = Group.1)) +
  #     ylab("") +
  #     xlab("Life history strategy") +
  #     ggtitle(title) +
  #     geom_line(color = "black") +
  #     geom_point(shape = 19) +
  #     facet_wrap(~variable, scales = "free") +
  #     theme_bw()
  #
  #   file_name <- paste0("figures/prelim_results/", climatic_variable, ".pdf")
  #   ggsave(file_name, height = 8, width = 14, units = "in")
  # }
}

# names(big_list) <- climatic_variables
# all_model_tables <- big_list
# all_aic_tables <- all_model_tables
# save(all_aic_tables, file = "all_aic_tables.Rsave")
# save(all_model_tables, file = "all_model_tables.Rsave")

# summ_tables <- do.call(rbind, lapply(climatic_variables, function(x) getSummTable(big_list, x)))
# write.csv(summ_tables, file = "tables/summ_table.csv")
# tmp <- do.call(rbind, lapply(all_aic_tables, function(x) do.call(rbind, x)))
# write.csv(tmp, file = "tables/model-weights.csv")
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # run  2 # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #

load("all_aic_tables.Rsave")
load("all_model_tables.Rsave")
phy <- read.tree("backbone_tree.tre")
phy$tip.label <- gsub("-.*", "", phy$tip.label)
phy <- drop.tip(phy, "Chorisporeae")
climatic_variables <- c("bio_1", "bio_4", "bio_5", "bio_6", "bio_12", "bio_14", "bio_15", "bio_ai")

# all_aic_tables$bio_1

# data organization
big_table <- do.call(rbind, lapply(all_aic_tables, function(x) do.call(rbind, x)))
big_table <- round(big_table, 2)
head(big_table)

tmp1 <- do.call(rbind, strsplit(rownames(big_table), "\\."))
tmp2 <- cbind(tmp1[,c(1,2)], do.call(rbind, strsplit(tmp1[,3], "_")))
tmp3 <- data.frame(variable = tmp2[,1], clade = tmp2[,2], model_class = tmp2[,3], model_type = tmp2[,4], AICcwt = big_table[,7])
total_support_table <- aggregate(tmp3$AICcwt, by = list(tmp3$model_class, tmp3$clade, tmp3$variable), sum)
# the values for our table are % character dependence withh some having mixtures
total_support_table$x <- total_support_table$x * c(1, 0, 0.5) 
plot_data <- aggregate(total_support_table$x, by = list(total_support_table$Group.2, total_support_table$Group.3), sum)

colnames(plot_data) <- c("id", "climate_variable", "percent_cd")
plot_data$climate_variable <- factor(plot_data$climate_variable, levels = c("bio_1", "bio_12", "bio_ai", "bio_4", "bio_15", "bio_5", "bio_6", "bio_14"))
plot_data <- plot_data[!plot_data$id == "Chorisporeae",]
# plotting
# ggplot(cd_support_table, aes(x = climate_variable, y =  id, fill = AICwt)) +
#   geom_tile()

a <- ggtree(phy) +
  geom_tiplab() +
  coord_cartesian(xlim = c(0, 160)) +
  ggtitle("a) Backbone Phylogeny")

b <- ggplot(plot_data, aes(x = climate_variable, y = id, fill = percent_cd)) +
  ggtitle("b) Support for character dependence") +
  geom_tile() + 
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme_tree2()

ab <- b %>% insert_left(a, width = 3)  
ggsave(filename = "figures/support_for_cd.pdf", plot = ab, height = 10, width = 15, units = "in")

library(dplyr)
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # individual plots # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
cols <- c("#8c510a", "#5ab4ac")
# variable <- climatic_variables[1]
a <- makePlot(climatic_variables[1], "a", "mean", TRUE)
b <- makePlot(climatic_variables[2], "b", "mean", TRUE)
c <- makePlot(climatic_variables[3], "c", "mean", TRUE)
d <- makePlot(climatic_variables[4], "d", "mean", TRUE)
e <- makePlot(climatic_variables[5], "e", "mean", TRUE)
f <- makePlot(climatic_variables[6], "f", "mean", TRUE)
g <- makePlot(climatic_variables[7], "g", "mean", TRUE)
h <- makePlot(climatic_variables[8], "h", "mean", TRUE)

ttest_table <- rbind(a, b, c, d, e, f, g, h)
write.csv(ttest_table, "tables/mean_ttest_table.csv")

a <- makePlot(climatic_variables[1], "a", "var", TRUE)
b <- makePlot(climatic_variables[2], "b", "var", TRUE)
c <- makePlot(climatic_variables[3], "c", "var", TRUE)
d <- makePlot(climatic_variables[4], "d", "var", TRUE)
e <- makePlot(climatic_variables[5], "e", "var", TRUE)
f <- makePlot(climatic_variables[6], "f", "var", TRUE)
g <- makePlot(climatic_variables[7], "g", "var", TRUE)
h <- makePlot(climatic_variables[8], "h", "var", TRUE)

ttest_table <- rbind(a, b, c, d, e, f, g, h)
write.csv(ttest_table, "tables/var_ttest_table.csv")

<<<<<<< HEAD
=======
final_plot <- grid.arrange(a, b, c, d, e, f, g, h, nrow=4)
ggsave("figures/mean-ttests.pdf", final_plot, height = 11, width = 8, units ="in")

a <- makePlot(climatic_variables[1], "a", "var")
b <- makePlot(climatic_variables[2], "b", "var")
c <- makePlot(climatic_variables[3], "c", "var")
d <- makePlot(climatic_variables[4], "d", "var")
e <- makePlot(climatic_variables[5], "e", "var")
f <- makePlot(climatic_variables[6], "f", "var")
g <- makePlot(climatic_variables[7], "g", "var")
h <- makePlot(climatic_variables[8], "h", "var")
>>>>>>> 5d09be531e8f7e605fce241453d05159017ccf10

final_plot <- grid.arrange(a, b, c, d, e, f, g, h, nrow=4)
ggsave("~/2022_life-history/figures/mean-ttests.pdf", final_plot, height = 10, width = 13, units = 
         "in")
ggsave("figures/var-ttests.pdf", final_plot, height = 11, width = 8, units ="in")

# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # individial analysis # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
climatic_variables
variable <- "bio_ai"

focal_list <- all_model_tables[[variable]]
focal_dat <- do.call(rbind, focal_list)
taxon_names <- unlist(lapply(focal_list, function(x) rownames(x)))
clade_names <- gsub("\\..*", "", rownames(focal_dat))
rownames(focal_dat) <- NULL
focal_dat$taxon <- taxon_names
focal_dat$clade <- clade_names
summ_data <- aggregate(focal_dat[,c(5,6)], by = list(focal_dat$clade, focal_dat$tip_state), mean)
summ_data[,3] <- convertVariable(variable, summ_data[,3])
summ_data <- summ_data[!summ_data$Group.1 == "Chorisporeae",] # remove Chorisporeae

diff_table <- aggregate(summ_data[,c(3,4)], by= list(summ_data$Group.1), diff) # perennial - annual
diff_table[match(range(diff_table[,2]), diff_table[,2]),]
colMeans(diff_table[,c(2,3)]) # perennial -  annual
diff_table[diff_table[,2] < 0,]
# 100 - 94
# 94 - 100

# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # other things (rates) # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #

big_table <- do.call(rbind, lapply(all_model_tables, function(x) do.call(rbind, x)))
aggregate(big_table$rates, by = list(big_table$tip_state), mean)
aggregate(big_table$rates, by = list(big_table$tip_state), sd)
range(big_table[big_table$tip_state == "annual", 1])
range(big_table[!big_table$tip_state == "annual", 1])


