
Codes used in: Long-term responses of life-history strategies to climatic variability in flowering plants (submitted)

Authors:
James D. Boyko, Eric R. Hagen, Jeremy M. Beaulieu, and Thais Vasconcelos




----
Description of folders: 
 
- **climate_data/** 

Values based on eight climatic variables for each occurrence point and summarized for each species in the dataset.

- **data_description/** 

Tables describing sampling fractions, number of tips and proportions of annuals and perennials for each analyzed clade.

- **datasets_final_for_hOUwie/** 

Datasets of climatic variables and life history strategies for each analyzed clade in the preferred input format for hOUwie.

- **tables/** 

Tables containing all final results from hOUwie analyses.

- **trait_dataset_post_curation/** 

Datasets of life history strategies post manual curation to include data beyond WCVP (source for included data is available in column "source" of each table)

- **trait_dataset_pre_curation/** 

Datasets of life history strategies as retrieved from the WCVP dataset.

- **trees_gbif_tips/** 

Phylogenetic trees with tips named according to the GBIF taxonomy backbone.

- **trees_simplified_tips/** 

Same as above, but without taxon authority. 

----
Scripts:

> 00_life_forms_for_manual_curation.R

Loads raw WCVP data and create a table for checking life forms.

> 01_taxizing.R

Creates reference table to match taxonomic backbones across all datasets (WCVP, GBIF, phylogenies)

> 02_visual_inspection.R

Script to plot tip states on the phylogenies to inspect distribution of life forms in all groups.

> 02.1_visual_inspection_post_curation.R

Same as above, but for trees that have passed manual curation (i.e. literature review of life forms).

> 03_send_names_to_gbif.R

Creates a occurrence dataset on GBIF for all species included in analyses following specific settings.

> 04_cleaning_gbif.R

Filters GBIF occurrence data for common inaccuracies in coordinates.

> 05_map_life_form.R

Creates maps of species richness and proportion of annual species as seen in Figure 1 of the main manuscript.

> 06_get_climate_vars.R

Creates dataset of means and infraspecific variation for each species and selected climatic variables imported from CHELSA.

> 07_organize_datasets_for_houwie.R

Organizes all datasets of traits and climatic variables so that they are in the preferred input format for hOUwie.

> 08_data-visualization.R

Creates plots matching phylogenetic trees, climatic variables, and character states per species for data visualization. 

> 09_run-houwie.R

Runs hOUwie models.

> 10-summarize-houwie.R

Summarizes results from hOUwie models.

> 11-backbone-tree.R

Creates a backbone tree where each tip is one of the analyzed clades based on the all inclusive tree of Smith and Brown (2018).

> 12-summarizeASR.R

Summarizes results from ancestral state reconstructions from hOUwie.

----
Other files:

> WCVP_life_form_description.txt

Description of life forms as presented by the WCVP website (https://powo.science.kew.org/about-wcvp#lifeforms).

----
Distribution data used in niche analyses comes from:

GBIF.org (29 April 2022) GBIF Occurrence Download https://doi.org/10.15468/dl.vt49kn

