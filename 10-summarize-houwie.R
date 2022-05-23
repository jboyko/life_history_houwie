require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(data.table)


setwd("~/2022_life-history/")

res_files <- dir("res_files/", full.names = TRUE)
group_names <- gsub("_bio.*", "", dir("res_files/"))

data_files <- dir("datasets_final_for_hOUwie/full_datasets/", full.names = TRUE)

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

all_model_fits <- lapply(group_names, quickSum)
names(all_model_fits) <- group_names
all_model_fits <- all_model_fits[!unlist(lapply(all_model_fits, is.null))]

model_tables <- lapply(all_model_fits, getModelTable)
lapply(model_tables, function(x) rownames(x)[x$AICwt > 1e-2])

all_model_avg_params <- lapply(all_model_fits, function(x) getModelAvgParams(x[1:8])[1:2])

model_set <- all_model_fits$Lupinus


odd_ball <- quickSum(group_names[4])

# Draba streptobrachia
# Draba brachystylis
odd_ball$CID_OUMA$data
require(phytools)
plotSimmap(OUwie:::correct_map_edges(odd_ball$CID_OUMA$simmaps[[5]]), fsize = 0.1)

pars <- odd_ball$CID_OUMA$p
mod <- odd_ball$CD_OUM
# hOUwie(phy = mod$phy, data = mod$data, rate.cat = 2, discrete_model = mod$index.disc, continuous_model = mod$index.cont, nSim = 25, p = pars)

getModelTable(odd_ball)

all_model_fits <- lapply(group_names, quickSum)
names(all_model_fits) <- group_names
all_model_fits <- all_model_fits[!unlist(lapply(all_model_fits, is.null))]

lapply(all_model_fits, getModelTable)

all_model_avg_params <- lapply(all_model_fits, function(x) getModelAvgParams(x)[1:2])


lapply(odd_ball, "[[", "simmaps")
all_maps <- lapply(odd_ball, "[[", "simmaps")
cd_maps <- all_maps[1:10]
cd_maps <- c(cd_maps[[1]], cd_maps[[2]], cd_maps[[3]], cd_maps[[4]], cd_maps[[5]],
             cd_maps[[6]], cd_maps[[7]], cd_maps[[8]], cd_maps[[9]], cd_maps[[10]])
cid_maps <- all_maps[11:18]
cid_maps <- c(cid_maps[[1]], cid_maps[[2]], cid_maps[[3]], cid_maps[[4]],
              cid_maps[[5]], cid_maps[[6]], cid_maps[[7]], cid_maps[[8]])

models <- c("BMV", "OUA", "OUV", "OUM", "OUMV", "OUMA", "OUVA", "OUMVA")

houwie_obj <- model_set[[1]]

new_maps <- lapply(model_set[1:8], "[[", "simmaps")
unlist(new_maps)


new_model_list <- mclapply(models, function(x) quickFunc(x), mc.cores = 8)


pars <- new_model_list$OUVA$p[c(1:11, 11)]
cont_model <- getOUParamStructure(model, 2, 2, TRUE)
res <- hOUwie.fixed(simmaps = cid_maps, data = odd_ball$CD_OUM$data, rate.cat = 2, discrete_model = odd_ball$CID_OUM$index.disc, continuous_model = "OUMVA", adaptive_sampling = FALSE, make_numeric = FALSE, diagn_msg = TRUE, p = pars)



group_name <- "Grewioideae"

dat_list <- organizeData(group_name, "bio_5", data_files, tree_files)
min(dat_list$phy$edge.length)

load("~/Desktop/weird-maps.Rsave")
new_maps <- lapply(simmaps, OUwie:::correct_map_edges)
pars <- odd_ball$CID_OUMA$p
mod <- odd_ball$CID_OUMA
map <- new_maps[[26]]
map_dat <- unlist(lapply(map$maps[map$edge[,2] <= Ntip(mod$phy)], function(x) as.numeric(names(x))[length(x)]))
map_sps <- mod$phy$tip.label[map$edge[,2][map$edge[,2] <= Ntip(mod$phy)]]
map_con <- mod$data[,3][match(map_sps, mod$data[,1])]
new_dat <- data.frame(sp = map_sps, reg = map_dat, x = map_con)

plotSimmap(new_maps[[26]], fsize = 0.1)

OUwie.fixed(phy = new_maps[[26]], data = new_dat, model = "OUMVA", simmap.tree = TRUE, 
            alpha = mod$solution.cont[1,], sigma.sq = mod$solution.cont[2,], theta = mod$solution.cont[3,])

adpt <- hOUwie(phy = mod$phy, data = mod$data, rate.cat = 2, discrete_model = mod$index.disc, continuous_model = mod$index.cont, nSim = 25, p = mod$p, adaptive_sampling = TRUE)

non_adpt <- hOUwie(phy = mod$phy, data = mod$data, rate.cat = 2, discrete_model = mod$index.disc, continuous_model = mod$index.cont, nSim = 25, p = mod$p, adaptive_sampling = FALSE)


adpt$expected_vals[names(adpt$expected_vals) == new_dat[,1][new_dat[,2] == 2]]
adpt$expected_vals[names(adpt$expected_vals) == new_dat[,1][new_dat[,2] == 1]]


mod$data


init_model_table <- getModelTable(odd_ball)

data.frame(TotlLik = init_model_table$lnLik, dTotlLik = init_model_table$lnLik - max(init_model_table$lnLik),
           DiscLik = init_model_table$DiscLik, dDiscLik = init_model_table$DiscLik - max(init_model_table$DiscLik),
           ContLik = init_model_table$ContLik, dContLik = init_model_table$ContLik - max(init_model_table$ContLik),
           row.names = rownames(init_model_table))

pars <- mod$p
pars[length(pars)-1] <- 5.5


good_maps <- new_maps[26:50]
for(i in 1:length(good_maps)){
  focal_map <- good_maps[[i]]$maps
  for(j in 1:length(focal_map)){
    names(focal_map[[j]])
  }
}


# pars[]
pars <- odd_ball$CID_OUMA$p[-8]
cid_oum_new <- hOUwie.fixed(simmaps = good_maps, data = mod$data, rate.cat = 2, discrete_model = odd_ball$CID_OUM$index.disc, continuous_model = odd_ball$CID_OUM$index.cont, ip = pars, adaptive_sampling = FALSE, make_numeric = FALSE, diagn_msg = TRUE)

pars <- odd_ball$CID_OUMA$p[c(1:6, 7, 8, 9, 9, 10, 11)]
cid_oumva_new <- hOUwie.fixed(simmaps = good_maps, data = mod$data, rate.cat = 2, discrete_model = odd_ball$CID_OUMVA$index.disc, continuous_model = odd_ball$CID_OUMVA$index.cont, ip = pars, adaptive_sampling = FALSE, make_numeric = FALSE, diagn_msg = TRUE)


mod$index.disc




