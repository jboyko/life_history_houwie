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

# Make sure the WCVP tables are in the same folder and load them
dist_sample <- read.csv("wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|")
names_sample <- read.csv("wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|")

# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")

# Now getting the list of species in your WCVP table
species_list <- unique(all_vars$taxon_name)

taxized_names_wcvp <- resolveGBIF(species_list) # This function adjust the names to the GBIF taxonomic backbone

# Make sure WCVP and GBIF communicate
reference_table <- data.frame(wcvp_name = species_list, gbif_name = taxized_names_wcvp) # you will need this table later
write.csv(reference_table, file="reference_table.csv", row.names = F) # saving table that you will need later

