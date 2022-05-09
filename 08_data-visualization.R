# imports
require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(geiger)
require(data.table)
require(reshape2)
require(ggplot2)
require(ggplotify)
require(gridExtra)
require(ggtree)
require(aplot)

# functions
organizeData <- function(clade_name, climate_variable, data_files, tree_files){
  focal_tree_file <- tree_files[grep(clade_name, tree_files)]
  focal_data_file <- data_files[grep(clade_name, data_files)]
  focal_data_file <- focal_data_file[grep(paste0(climate_variable, "_"), focal_data_file)]
  phy <- read.tree(focal_tree_file)
  dat <- read.csv(focal_data_file)
  dat$species <- gsub(" ", "_", dat$species)
  dat <- dat[match(phy$tip.label, dat$species),]
  cat("\n", "Removed", length(which(dat$life_form == "no_life_form_on_database")), "species of", length(dat$species), "because they didn't have life history data.\n")
  dat <- dat[!dat$life_form == "no_life_form_on_database",]
  plot_data <- data.frame(id = dat$species, value = dat[,"mean"], life_form = as.factor(dat[,"life_form"]))
  cat("\n", "Removed", length(which(apply(plot_data, 1, function(x) !any(is.na(x))))), "species of", length(dat$species), "because they didn't have climate data.\n")
  plot_data <- plot_data[apply(plot_data, 1, function(x) !any(is.na(x))),]
  pruned_phy <- keep.tip(phy, phy$tip.label[match(plot_data$id, phy$tip.label)])
  return(list(dat=plot_data[,c(1,3,2)], phy = pruned_phy))
}


makePlot <- function(clade_name, climatic_variable, data_files, tree_files){
  dat_list <- organizeData(clade_name, climatic_variable, data_files, tree_files)
  phy <- dat_list$phy
  dat <- dat_list$dat
  dat[,3] <- log(dat[,3])
  a <- ggtree(phy) +
    ggtitle(clade_name)
  b <- ggplot(dat, aes(x = id, y = value, color = life_form, fill = life_form)) + 
    geom_col(width = .000001) + 
    scale_color_manual(values = cols_per_state[1:2]) +
    scale_fill_manual(values = cols_per_state[1:2]) +
    coord_flip() + 
    theme_tree2() + 
    ggtitle(climatic_variable)
  
  ab <- b %>% insert_left(a, width = 3)  
  # return(ab)
  file_name <- paste0("figures/preim_raw_data/", clade_name, "_", climatic_variable, ".pdf")
  ggsave(file_name, ab, height = 8, width = 12, units = "in")
}

runTtest <- function(clade_name, climatic_variable, data_files, tree_files){
  dat_list <- organizeData(clade_name, climatic_variable, data_files, tree_files)
  dat <- dat_list$dat[,c(2,3)]
  # dat[,2] <- log(dat[,2])
  model <- t.test(value ~ life_form, data = dat)
  
  out <- data.frame(clade_name, climatic_variable, model$statistic, model$p.value, model$estimate[1], model$estimate[2], model$parameter)
  colnames(out) <- c("clade", "climate", "t-stat", "p-value", "annual", "perennial", "df")
  rownames(out) <- NULL
  return(out)
}

# working directory
setwd("~/2022_life-history/")

data_files <- dir("datasets_final_for_hOUwie/full_datasets/", full.names = TRUE)
group_names <- unique(unlist(lapply(strsplit(dir("datasets_final_for_hOUwie/full_datasets/"), "-"), function(x) x[[1]])))

tree_files <- dir("trees_simplified_tips/", full.names = TRUE)
clade_name <- group_names[1]

# sapply(group_names, function(x) makePlot(x, "bio_15", data_files, tree_files))

climatic_variables <- c(paste0("bio_", 1:19), "bio_ai", "bio_et0")

big_df <- do.call(rbind, lapply(climatic_variables, function(y) do.call(rbind, lapply(group_names, function(x) runTtest(x, y, data_files, tree_files)))))

aggregate(big_df[,c(4,5,6)], by = list(big_df$climate), median)
aggregate(big_df[,c(4,5,6)], by = list(big_df$clade), mean)

write.csv(big_df, file = "prelim/all_t-tests.csv")

# pdf("figures/preim_raw_data/bio_15.pdf", onefile = TRUE)
for(i in 1:length(group_names)){
  makePlot(group_names[i], "bio_5", data_files, tree_files)
}
# dev.off()

organizeData(group_names[i], "bio_5", data_files, tree_files)
hist(log(dat_list$dat[,3]))


cols_per_state <- c("#a6611a", "#018571", "#dfc27d", "#80cdc1")


phy <- dat_list$phy
dat <- dat_list$dat
dat[,3] <- log(dat[,3])

a <- ggtree(phy) +
  ggtitle(clade_name)
b <- ggplot(dat, aes(x = id, y = value, color = life_form, fill = life_form)) + 
  geom_col(width = .000001) + 
  scale_color_manual(values = cols_per_state[1:2]) +
  scale_fill_manual(values = cols_per_state[1:2]) +
  coord_flip() + 
  theme_tree2() + 
  ggtitle(climatic_variable)

ab <- b %>% insert_left(a, width = 3)  
ab

file_name <- paste0("figures/preim_raw_data/", clade_name, "_", climatic_variable, ".pdf")
ggsave(file_name, ab, height = 8, width = 12, units = "in")

out_1 <- hOUwie(phy, dat, 1, "ER", "OUM", nSim = 25, diagn_msg = TRUE, ub_continuous_model = c(1e10, 1e10, 1e10))
out_2 <- hOUwie(phy, dat, 1, "ER", "OU1", nSim = 25, diagn_msg = TRUE, ub_continuous_model = c(1e10, 1e10, 1e10))
out_3 <- hOUwie(phy, dat, 1, "ER", "BM1", nSim = 25, diagn_msg = TRUE, ub_continuous_model = c(1e10, 1e10, 1e10))
out_4 <- hOUwie(phy, dat, 2, "ER", "OUM", nSim = 25, null.model = TRUE, diagn_msg = TRUE)


OUwie(phy, dat, "BM1")

Q <- matrix(c(-1,1,1,-1)/100, 2, 2)
map <- makeSimmap(phy, dat[,c(1,2)], Q, 1)[[1]]


OUwie(map, dat, "OUM", simmap.tree = TRUE, starting.vals = c(1, 50), ub = c(1e10))

plot(map)


