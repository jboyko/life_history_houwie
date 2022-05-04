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
organizeData <- function(clade_name, data_files, tree_files){
  focal_tree_file <- tree_files[grep(clade_name, tree_files)]
  focal_data_file <- data_files[grep(clade_name, data_files)]
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

# working directory
setwd("2022_life-history/")
# setwd("~/Desktop/WCVP_special_issue/James_perennial_annual/life_history_houwie")

data_files <- dir("datasets_final_for_hOUwie/full_datasets", full.names = TRUE)
group_names <- unlist(lapply(strsplit(dir("datasets_final_for_hOUwie/full_datasets/"), "-"), function(x) x[[1]]))
tree_files <- dir("trees_simplified_tips/", full.names = TRUE)
clade_name <- group_names[11]

dat_list <- organizeData(group_names[11], data_files, tree_files)

cols_per_state <- c("#a6611a", "#018571", "#dfc27d", "#80cdc1")

phy <- dat_list$phy
dat <- dat_list$dat

a <- ggtree(phy)
b <- ggplot(dat, aes(x = id, y = value, color = life_form)) + 
  geom_col(width = .000001) + 
  scale_color_manual(values = cols_per_state[1:2]) +
  coord_flip() + 
  theme_tree2() + 
  ylab("Bio15") +
  ggtitle("")

ab <- b %>% insert_left(a, width = 3)  
ab


out_1 <- hOUwie(phy, dat, 1, "ER", "OUM", nSim = 25, diagn_msg = TRUE, ub_continuous_model = c(1e10, 1e10, 1e10))
out_2 <- hOUwie(phy, dat, 1, "ER", "OU1", nSim = 25, diagn_msg = TRUE, ub_continuous_model = c(1e10, 1e10, 1e10))
out_3 <- hOUwie(phy, dat, 1, "ER", "BM1", nSim = 25, diagn_msg = TRUE, ub_continuous_model = c(1e10, 1e10, 1e10))
out_4 <- hOUwie(phy, dat, 2, "ER", "OUM", nSim = 25, null.model = TRUE, diagn_msg = TRUE)


OUwie(phy, dat, "BM1")

Q <- matrix(c(-1,1,1,-1)/100, 2, 2)
map <- makeSimmap(phy, dat[,c(1,2)], Q, 1)[[1]]


OUwie(map, dat, "OUM", simmap.tree = TRUE, starting.vals = c(1, 50), ub = c(1e10))

plot(map)


