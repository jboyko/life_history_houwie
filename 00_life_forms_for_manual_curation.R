#----------------
dist_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")
all_life_forms <- all_vars$lifeform_description
all_life_forms <- as.data.frame(table(all_life_forms))
all_life_forms <- all_life_forms[-1,]
all_life_forms$scoring1 <- NA
write.csv(all_life_forms, file=paste0(Sys.Date(),"_life_form.csv"), row.names=F)
