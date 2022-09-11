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

makePlot <- function(variable, letter, mu){
  climate_variable <- getVariableName(variable)
  
  focal_list <- all_model_tables[[variable]]
  focal_dat <- do.call(rbind, focal_list)
  taxon_names <- unlist(lapply(focal_list, function(x) rownames(x)))
  clade_names <- gsub("\\..*", "", rownames(focal_dat))
  rownames(focal_dat) <- NULL
  focal_dat$taxon <- taxon_names
  focal_dat$clade <- clade_names
  # head(focal_dat)
  summ_data <- aggregate(focal_dat[,c(5,6)], by = list(focal_dat$clade, focal_dat$tip_state), mean)
  if(mu == "mean" | mu == "both"){
    n_annual <- length(which(summ_data[summ_data$Group.2 == "annual", 3] - summ_data[summ_data$Group.2 == "perennial", 3] > 0))
    ttest_res <- runTtest(phy, summ_data, 3)
    if(ttest_res$dbar < 0){
      n_annual <- 33 - n_annual
    }
    summ_data[,3] <- convertVariable(variable, summ_data[,3])
    a <- ggplot(summ_data, aes(x = Group.2, y = expected_mean, group = Group.1, color = Group.2)) +
      xlab("") +
      ylab("Expected mean") +
      geom_line(color = "light grey") +
      geom_point(shape = 19, color = "light grey") +
      ggtitle(paste0(letter, ") ", climate_variable), subtitle = paste0("p = ", round(ttest_res$P.dbar, 3), " (", n_annual, " out of 33 clades)" )) +
      theme_bw() +
      stat_summary(fun=mean, geom="point",aes(group=1, size = 2), color = cols) +  
      stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group=1), width = 0.15, color = "black") +
      # labs(caption = paste0("p = ", round(t_sigma$P.dbar, 3))) +
      theme(legend.position="none", text = element_text(size = 15), axis.text.y = element_text(size = 10))
    out <- a
  }
  if(mu == "var" | mu == "both"){
    n_annual <- length(which(summ_data[summ_data$Group.2 == "annual", 4] - summ_data[summ_data$Group.2 == "perennial", 4] > 0))
    ttest_res <- runTtest(phy, summ_data, 4)
    if(ttest_res$dbar < 0){
      n_annual <- 33 - n_annual
    }
    b <- ggplot(summ_data, aes(x = Group.2, y = expected_var, group = Group.1, color = Group.2)) +
      xlab("") +
      ylab("Expected mean") +
      ggtitle(paste0(letter, ") ", climate_variable), subtitle = paste0("p = ", round(ttest_res$P.dbar, 3), " (", n_annual, " out of 33 clades)" )) +
      geom_line(color = "light grey") +
      geom_point(shape = 19, color = "light grey") +
      theme_bw() +
      stat_summary(fun=mean, geom="point",aes(group=1, size = 2), color = cols) +  
      stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group=1), width = 0.15, color = "black") +
      # labs(caption = paste0("p = ", round(t_sigma$P.dbar, 3))) +
      theme(legend.position="none", text = element_text(size = 15), axis.text.y = element_text(size = 10))
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
# # # # # run  1 # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #

climatic_variables <- c("bio_1", "bio_4", "bio_5", "bio_6", "bio_12", "bio_14", "bio_15", "bio_ai")
# climatic_variables <- c("bio_6")
big_list <- list()
print <- FALSE

# for(i in 1:length(climatic_variables)){
#   climatic_variable <- climatic_variables[i]
#   print("loading")
#   mod_avg_param_list <- lapply(group_names, function(x) get_mod_avg_param_table(climatic_variable, x, type="AICc"))
#   # test <- load_file(climatic_variable, "Balsamiaceae")
#   # mod_avg_param_list <- lapply(group_names, function(x) get_mod_table(climatic_variable, x))
#   names(mod_avg_param_list) <- group_names
#   big_list[[i]] <- mod_avg_param_list
#   # get_mod_table(climatic_variable, "CES")
#   # test <- mod_avg_param_list[[1]]
#   # mod_avg_param_list <- big_list[[i]]
#   # aggregate(test[,1:4], by = list(test$tip_state), mean)
#   # if(print){
#   #   plot_data <- do.call(rbind, lapply(mod_avg_param_list, function(x) aggregate(x[,1:4], by = list(x$tip_state), mean)))
#   #   data_se <- do.call(rbind, lapply(mod_avg_param_list, function(x) aggregate(x[,1:4], by = list(x$tip_state), get_stderr)))
#   #   plot_data$clade <- gsub("\\..*", "", rownames(plot_data))
#   #   title <- getVariableName(climatic_variable)
#   #   
#   #   # discrete rate
#   #   plot_data_1 <- plot_data[,c(1,2,6)]
#   #   ranges <- c(mean(plot_data_1[,2]) - (sd(plot_data_1[,2])), mean(plot_data_1[,2]) + (sd(plot_data_1[,2])))
#   #   p1 <- ggplot(melt(plot_data_1), aes(x = Group.1, y = value, group = clade, color = Group.1)) +
#   #     ylab("") +
#   #     xlab("Life history strategy") +
#   #     ggtitle("Discrete rate") +
#   #     geom_line(color = "black") +
#   #     geom_point(shape = 19) +
#   #     coord_cartesian(ylim=ranges) +
#   #     theme_bw()
#   #   
#   #   # alpha
#   #   plot_data_1 <- plot_data[,c(1,3,6)]
#   #   ranges <- c(mean(plot_data_1[,2]) - (sd(plot_data_1[,2])), mean(plot_data_1[,2]) + (sd(plot_data_1[,2])))
#   #   p2 <- ggplot(melt(plot_data_1), aes(x = Group.1, y = value, group = clade, color = Group.1)) +
#   #     ylab("") +
#   #     xlab("Life history strategy") +
#   #     ggtitle("Alpha") +
#   #     geom_line(color = "black") +
#   #     geom_point(shape = 19) +
#   #     coord_cartesian(ylim=ranges) +
#   #     theme_bw()
#   #   
#   #   # sigma.sq
#   #   plot_data_1 <- plot_data[,c(1,4,6)]
#   #   ranges <- c(mean(plot_data_1[,2]) - (sd(plot_data_1[,2])), mean(plot_data_1[,2]) + (sd(plot_data_1[,2])))
#   #   p3 <- ggplot(melt(plot_data_1), aes(x = Group.1, y = value, group = clade, color = Group.1)) +
#   #     ylab("") +
#   #     xlab("Life history strategy") +
#   #     ggtitle("Sigma Squared") +
#   #     geom_line(color = "black") +
#   #     geom_point(shape = 19) +
#   #     coord_cartesian(ylim=ranges) +
#   #     theme_bw()
#   #   
#   #   # sigma.sq
#   #   plot_data_1 <- plot_data[,c(1,5,6)]
#   #   ranges <- c(mean(plot_data_1[,2]) - (sd(plot_data_1[,2])), mean(plot_data_1[,2]) + (sd(plot_data_1[,2])))
#   #   p4 <- ggplot(melt(plot_data_1), aes(x = Group.1, y = value, group = clade, color = Group.1)) +
#   #     ylab("") +
#   #     xlab("Life history strategy") +
#   #     ggtitle("Theta") +
#   #     geom_line(color = "black") +
#   #     geom_point(shape = 19) +
#   #     coord_cartesian(ylim=ranges) +
#   #     theme_bw()
#   #   
#   #   # plot_data <- melt(plot_data, by = list("waiting_times", "alpha", "sigma.sq", "theta"))
#   # 
#   #   title <- getVariableName(climatic_variable)
#   #   print("plotting")
#   #   ggplot(plot_data, aes(x = Group.1, y = value, group = clade, color = Group.1)) +
#   #     ylab("") +
#   #     xlab("Life history strategy") +
#   #     ggtitle(title) +
#   #     geom_line(color = "black") +
#   #     geom_point(shape = 19) +
#   #     facet_wrap(~variable, scales = "free") +
#   #     theme_bw()
#   # 
#   #   file_name <- paste0("figures/prelim_results/", climatic_variable, ".pdf")
#   #   ggsave(file_name, height = 8, width = 14, units = "in")
#   # }
# }

# names(big_list) <- climatic_variables
# all_model_tables <- big_list
# save(all_model_tables, file = "all_model_tables.Rsave")

# summ_tables <- do.call(rbind, lapply(climatic_variables, function(x) getSummTable(big_list, x)))
# write.csv(summ_tables, file = "tables/summ_table.csv")

# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # run  2 # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #

load("all_model_tables.Rsave")
phy <- read.tree("backbone_tree.tre")
phy$tip.label <- gsub("-.*", "", phy$tip.label)
climatic_variables <- c("bio_1", "bio_4", "bio_5", "bio_6", "bio_12", "bio_14", "bio_15", "bio_ai")

for(i in 1:length(all_model_tables)){
  focal_tbl <- do.call(rbind, all_model_tables[[i]])
  write.csv(focal_tbl, file = paste0("tables/model_tables/model_table-", names(all_model_tables)[i], ".csv"))
}

# data organization
big_table <- do.call(cbind, lapply(lapply(all_model_tables, function(x) do.call(rbind, x)), function(y) round(y, 2)))
aicwt_table <- melt(cbind(rownames(big_table), big_table[,grep("AICwt", colnames(big_table))]))
tmp1 <- do.call(rbind, strsplit(aicwt_table[,1], "\\."))
tmp2 <- cbind(tmp1[,1], do.call(rbind, strsplit(tmp1[,2], "_")))
aicwt_table[,2] <- gsub("\\..*", "", aicwt_table[,2])
tmp3 <- data.frame(id = tmp2[,1], model_class = tmp2[,2], model_type = tmp2[,3], aicwt_table[,c(2,3)])
cd_support_table <- aggregate(tmp3$value, by = list(tmp3$id, tmp3$variable, tmp3$model_class), sum)
cd_support_table <- cd_support_table[cd_support_table$Group.3 == "CD",]
colnames(cd_support_table) <- c("id", "climate_variable", "model_class", "AICwt")
cd_support_table$climate_variable <- factor(cd_support_table$climate_variable, levels = c("bio_1", "bio_12", "bio_ai", "bio_4", "bio_15", "bio_5", "bio_6", "bio_14"))

# plotting
# ggplot(cd_support_table, aes(x = climate_variable, y =  id, fill = AICwt)) +
#   geom_tile()

a <- ggtree(phy) +
  geom_tiplab() +
  coord_cartesian(xlim = c(0, 160)) +
  ggtitle("a) Backbone Phylogeny")

b <- ggplot(cd_support_table, aes(x = climate_variable, y = id, fill = AICwt)) +
  ggtitle("b) Support for character dependence") +
  geom_tile() + 
  scale_fill_viridis_c() +
  theme_tree2()

ab <- b %>% insert_left(a, width = 3)  
ggsave(filename = "figures/support_for_cd.pdf", plot = ab, height = 10, width = 15, units = "in")


# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # individual plots # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
cols <- c("#8c510a", "#5ab4ac")
# variable <- climatic_variables[1]
a <- makePlot(climatic_variables[1], "a", "var")
b <- makePlot(climatic_variables[2], "b", "var")
c <- makePlot(climatic_variables[3], "c", "var")
d <- makePlot(climatic_variables[4], "d", "var")
e <- makePlot(climatic_variables[5], "e", "var")
f <- makePlot(climatic_variables[6], "f", "var")
g <- makePlot(climatic_variables[7], "g", "var")
h <- makePlot(climatic_variables[8], "h", "var")

final_plot <- grid.arrange(a, b, c, d, e, f, g, h, nrow=4)
ggsave("~/2022_life-history/figures/vars-ttests.pdf", final_plot, height = 10, width = 13, units = "in")

# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # some tip rate stuff  # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
load("all_model_avg_res.Rsave")
big_list <- all_model_avg_res

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

phy <- read.tree("backbone_tree.tre")
phy$tip.label <- gsub("-.*", "", phy$tip.label)
cols <- c("#8c510a", "#5ab4ac")


plot_list <- list()
for(i in 1:8){
  print(i)
  focal_list <- big_list[[i]]
  tmp_table <- do.call(rbind, lapply(focal_list, function(x) aggregate(x[,1:4], by = list(x$tip_state), mean)))
  tmp_table <- cbind(clade=gsub("\\..*", "", rownames(tmp_table)), tmp_table)
  # problems <- c("Gesneriaceae", "Balsamiaceae")
  # tmp_table <- tmp_table[-grep(problems[1], rownames(tmp_table)),]
  # tmp_table <- tmp_table[-grep(problems[2], rownames(tmp_table)),]
  rownames(tmp_table) <- NULL
  names(big_list)[i]
  diff_table <- aggregate(tmp_table[,3:6], by = list(tmp_table$clade), function(x) diff(exp(x)))
  diff_table[order(diff_table[,5]),]
  plot_data <- melt(tmp_table, by = list("waiting_times", "alpha", "sigma.sq", "theta"))
  
  title <- getVariableName(names(big_list)[i])
  
  # t_trates <- runTtest(phy, plot_data, 'waiting_times')
  # p_trates <- ggplot(subset(plot_data, plot_data$variable == 'waiting_times'), aes(x = Group.1, y = value, group = clade, color = Group.1)) +
  #   ylab("") +
  #   xlab("") +
  #   ggtitle("a) Transition rate") +
  #   geom_line(color = "light grey") +
  #   geom_point(shape = 19, color = "light grey") +
  #   theme_bw() +
  #   stat_summary(fun=mean,geom="point",aes(group=1, size = 2), color = cols) +  
  #   stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group=1), width = 0.15, color = "black") +
  #   labs(caption = paste0("p=",round(t_trates$P.dbar, 3))) + 
  #   theme(legend.position="none")
  
  t_sigma <- runTtest(phy, plot_data, 'sigma.sq')
  p_sigma <- ggplot(subset(plot_data, plot_data$variable == 'sigma.sq'), 
                    aes(x = Group.1, y = value, group = clade, color = Group.1)) +
    ylab(title) +
    xlab("") +
    ggtitle("", subtitle = paste0("p = ", round(t_sigma$P.dbar, 3))) +
    geom_line(color = "light grey") +
    geom_point(shape = 19, color = "light grey") +
    theme_bw() +
    stat_summary(fun=mean,geom="point",aes(group=1, size = 2), color = cols) +  
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group=1), width = 0.15, color = "black") +
    # labs(caption = paste0("p = ", round(t_sigma$P.dbar, 3))) +
    theme(legend.position="none", text = element_text(size = 15), axis.text.y = element_text(size = 10))
  
  t_theta <- runTtest(phy, plot_data, 'theta')
  p_theta <- ggplot(subset(plot_data, plot_data$variable == 'theta'), 
                    aes(x = Group.1, y = value, group = clade, color = Group.1)) +
    ylab("") +
    xlab("") +
    ggtitle("", subtitle = paste0("p = ", round(t_theta$P.dbar, 3))) +
    geom_line(color = "light grey") +
    geom_point(shape = 19, color = "light grey") +
    stat_summary(fun=mean,geom="point",aes(group=1, size = 2), color = cols) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group=1), width = 0.15, color = "black") +
    theme_bw() +
    # labs(caption = paste0("p = ", round(t_theta$P.dbar, 3))) +
    theme(legend.position="none", text = element_text(size = 15), axis.text.y = element_text(size = 10))
  
  
  plot_list[[i]] <- grid.arrange(p_sigma, p_theta, nrow=1)
}

poop1 <- grid.arrange(plot_list[[1]], plot_list[[5]], plot_list[[8]], ncol = 1)

poop2 <- grid.arrange(plot_list[[2]], plot_list[[7]], ncol = 1)

poop3 <- grid.arrange(plot_list[[3]], plot_list[[4]], plot_list[[6]], ncol = 1)


ggsave(filename = "figures/ttest-plots-mean.pdf",  plot = poop1, height = 15, width = 10, units = "in")
ggsave(filename = "figures/ttest-plots-vars.pdf",  plot = poop2, height = 10, width = 10, units = "in")
ggsave(filename = "figures/ttest-plots-extr.pdf",  plot = poop3, height = 15, width = 10, units = "in")


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

# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # stuff that summarizes clade specific stuff  # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #



# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # stuff that summarizes other stuff  # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #


phy <- read.tree("backbone_tree.tre")
phy$tip.label <- gsub("-.*", "", phy$tip.label)
climatic_variable <- "bio_1"


summ_tables <- do.call(rbind, lapply(climatic_variables, function(x) getSummTable(big_list, x)))
write.csv(summ_tables, file = "tables/summ_table.csv")

getSummTable(big_list, "bio_ai")









