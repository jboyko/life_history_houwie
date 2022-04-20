# Taxizing data
library(taxize)

#' Taxize GBIF
#' @param name A character vector with species names.
resolveGBIF <- function(name) {
  gnr_resolve_x <- function(x) {
    sources <- taxize::gnr_datasources()
    tmp.name <- suppressWarnings(taxize::gnr_resolve(names=x, data_source_ids=sources$id[sources$title == "GBIF Backbone Taxonomy"], best_match_only=TRUE)$matched_name)
    if(is.null(tmp.name)) {
      tmp.name <- paste0("UNMATCHED_",x)
    }
    return(tmp.name)
  }
  new.names <- pbapply::pbsapply(name, gnr_resolve_x)
  return(as.character(new.names))
}

taxize.batch <- function(species_names, file_name="") {
  taxized_names_wcvp <- resolveGBIF(species_names) 
  reference_table <- data.frame(wcvp_name = species_names, gbif_name = taxized_names_wcvp) 
  write.csv(reference_table, file=file_name, row.names = F)  
}


# Make sure the WCVP tables are in the same folder and load them
#dist_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
#names_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

dist_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")

# Now getting the list of species in your WCVP table
species_list <- unique(all_vars$taxon_name)
species_list <- subset(species_list, species_list!="")

species_list_1a <- species_list[1:30000]
species_list_1b <- species_list[30001:60000]
species_list_1c <- species_list[60001:100000]

species_list_2 <- species_list[100001:200000]
species_list_3 <- species_list[200001:300000]
species_list_4 <- species_list[300001:length(species_list)]

try(taxize.batch(species_list_1a, file_name="reference_table_1a.csv"))
try(taxize.batch(species_list_1b, file_name="reference_table_1b.csv"))
try(taxize.batch(species_list_1c, file_name="reference_table_1c.csv"))

try(taxize.batch(species_list_2, file_name="reference_table_2.csv"))
try(taxize.batch(species_list_3, file_name="reference_table_3.csv"))
try(taxize.batch(species_list_4, file_name="reference_table_4.csv"))


