# Setting table for manual curation of life forms:

setwd("~/Desktop/WCVP_special_issue/James_perennial_annual/life_history_houwie")
#rm(list=ls())

#----------------
# Load Kew data:
dist_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")
all_life_forms <- all_vars$lifeform_description

all_life_forms <- as.data.frame(table(all_life_forms))
all_life_forms <- all_life_forms[-1,]
all_life_forms$scoring1 <- NA
write.csv(all_life_forms, file=paste0(Sys.Date(),"_life_form.csv"), row.names=F)


#all_climbers <- subset(all_vars, grepl(paste0(c("Climb","Liana"), collapse="|"), all_vars$lifeform_description))
#write.csv(all_climbers, file="all_climbers.csv", row.names=F)
#length(unique(all_climbers$taxon_name))

#all_spp_climbers  <- unique(all_climbers$taxon_name)
#all_climbing_habit <- unique(all_climbers$lifeform_description)
#result <- data.frame(life_form=all_climbing_habit, n_species=NA)
#for(i in 1:length(all_climbing_habit)) {
#  all_climbers_tmp <- subset(all_climbers, all_climbers$lifeform_description==all_climbing_habit[i])
#  n_tmp <- unique(all_climbers_tmp$taxon_name)
#  result[i,2] <- length(n_tmp)
#}
#write.csv(result, file="n_species_climbers.csv", row.names = F)

#all_life_forms <- as.data.frame(table(all_climbers$lifeform_description))
#all_life_forms <- all_life_forms[-1,]
#all_life_forms$scoring1 <- NA
#write.csv(all_life_forms, file=paste0(Sys.Date(),"_life_form.csv"), row.names=F)