require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(data.table)


setwd("~/2022_life-history/")

res_files <- dir("res_files/", full.names = TRUE)

data_files <- dir("datasets_final_for_hOUwie/full_datasets/", full.names = TRUE)
group_names <- unique(unlist(lapply(strsplit(dir("datasets_final_for_hOUwie/full_datasets/"), "-"), function(x) x[[1]])))

tree_files <- dir("trees_simplified_tips/", full.names = TRUE)


quickSum <- function(group_name){
  load(res_files[grep(group_name, res_files)])
  if(all(unlist(lapply(model_set_res, class)) == "houwie")){
    out <- model_set_res
  }else{
    out <- NULL
  }
  return(out)
}
# out <- c()
# for(i in 1:length(group_names)){
#   group_name <- group_names[i]
#   load(res_files[grep(group_name, res_files)])
#   if(!all(unlist(lapply(model_set_res, class)) == "houwie")){
#     out <- rbind(out, c(group_name, rep(NA, 7)))
#   }else{
#     out <- rbind(out, )
#   }
# }

all_model_fits <- lapply(group_names, quickSum)
names(all_model_fits) <- group_names
all_model_fits <- all_model_fits[!unlist(lapply(all_model_fits, is.null))]

all_model_avg_params <- lapply(all_model_fits, function(x) getModelAvgParams(x)[1:2])




group_name <- "Grewioideae"

dat_list <- organizeData(group_name, "bio_5", data_files, tree_files)
min(dat_list$phy$edge.length)
