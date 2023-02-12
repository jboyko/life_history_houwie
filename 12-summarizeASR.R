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
get_mod_avg_recon <- function(climatic_variable, clade_name){
  focal_files <- dir(res_folders[grep(paste0(climatic_variable, "$"), res_folders)], full.names = TRUE)
  focal_file <- focal_files[grep(clade_name, focal_files)]
  load(focal_file)
  mod_avg_param_table <- getModelTable(complete_list, "AICc")
  complete_list <- complete_list[!is.na(names(complete_list))] # remove failed
  root_p_list <- lapply(complete_list, getRootFromModel)
  cd_index <- unlist(lapply(root_p_list, function(x) length(x) == 2))
  cd_roots <- do.call(rbind, root_p_list[cd_index])
  cid_roots <- do.call(rbind, root_p_list[!cd_index])[,c(1,2)] + do.call(rbind, root_p_list[!cd_index])[,c(3,4)]
  root_probs <- rbind(cd_roots, cid_roots)
  weighted_root_state <- colSums(root_probs * mod_avg_param_table[,7])
  weighted_root_state <- weighted_root_state/sum(weighted_root_state)
  return(weighted_root_state)
}

getRootStateFromMap <- function(map){
  return(names(map$maps[[which.min(map$edge[,1])]])[1])
}

get_complete_list <- function(climatic_variable, clade_name){
  focal_files <- dir(res_folders[grep(paste0(climatic_variable, "$"), res_folders)], full.names = TRUE)
  focal_file <- focal_files[grep(clade_name, focal_files)]
  load(focal_file)
  return(complete_list)
}

getRootFromModel <- function(model){
  root_matrix <- matrix(0, length(model$simmaps), dim(model$continuous_model)[2])
  colnames(root_matrix) <- colnames(model$continuous_model)
  root_states <- as.numeric(unlist(lapply(model$simmaps, getRootStateFromMap)))
  joint_probs <- model$all_cont_liks + model$all_disc_liks
  joint_weights <- exp(joint_probs - max(joint_probs))/sum(exp(joint_probs - max(joint_probs)))
  for(i in 1:length(root_states)){
    root_matrix[i,root_states[i]] <- 1 * joint_weights[i]
  }
  root_probs <- colSums(root_matrix)/sum(root_matrix)
  return(root_probs)
}

# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # setup # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #

setwd("~/2022_life-history/")

res_folders <- dir("res_files/compiled_models/", full.names = TRUE)
group_names <- unique(unlist(lapply(strsplit(dir("datasets_final_for_hOUwie/full_datasets/"), "-"), function(x) x[[1]])))

data_files <- dir("datasets_final_for_hOUwie/full_datasets/", full.names = TRUE)

tree_files <- dir("trees_simplified_tips/", full.names = TRUE)

climatic_variables <- c("bio_1", "bio_4", "bio_5", "bio_6", "bio_12", "bio_14", "bio_15", "bio_ai")

# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # # run  # # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # #

asr_list <- list()
complete_list <- list()
for(i in 1:length(group_names)){
  print(group_names[i])
  tmp_list <- list()
  for(j in 1:length(climatic_variables)){
    tmp_list[[j]] <- get_mod_avg_recon(climatic_variables[j], group_names[i])
  }
  asr_list[[i]] <- do.call(rbind, tmp_list)
  rownames(asr_list[[i]]) <- climatic_variables
}

names(asr_list) <- group_names

asr_table <- data.frame(clade = rep(group_names, each = length(climatic_variables)), 
                        variable = rep(climatic_variables, length(group_names)), 
                        do.call(rbind, asr_list), row.names = NULL)
write.csv(asr_table, file = "tables/asr_table.csv")

all_asrs <- do.call(rbind, lapply(asr_list, function(x) colMeans(x)/sum(colMeans(x))))
p_annual <- data.frame(clade = rownames(all_asrs), p_ann = all_asrs[,1], p_min = do.call(rbind, lapply(asr_list, function(x) range(x[,1])))[,1], p_max = do.call(rbind, lapply(asr_list, function(x) range(x[,1])))[,2])
p_annual <- p_annual[!p_annual$clade == "Chorisporeae",]

p <- ggplot(p_annual, aes(x = reorder(clade, p_ann), y = p_ann)) +
  geom_point() +
  geom_errorbar(aes(ymin=p_min, ymax=p_max), width=.2) +
  xlab("Clade") +
  ylab("Probability of annual root state") +
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip(ylim=c(0,1)) +
  theme_bw()

ggsave("figures/asr-plot.pdf", p, height = 10, width = 12, units = "in")


# asr_list$Apioideae
bio_5 <- get_complete_list("bio_5", "Apioideae")
bio_12 <- get_complete_list("bio_12", "Apioideae")

getModelTable(bio_5)
getModelTable(bio_12)
bio_12$CD_OUM
bio_5$CD_OUM
