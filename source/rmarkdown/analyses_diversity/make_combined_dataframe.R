# Originally based on Silke's code from issue #152
# comment https://github.com/slambrechts/INBO_eDNA_metabarcoding_BODEM/issues/152#issuecomment-2298667395

# Make combined dataframe for MBAG

# Use not-rarefied dataframes

library(dplyr)
library(phyloseq)
library(tidytacos)
library(tibble)

# Current workdir
print("Current working directory:")
print(getwd())

# Source R helper functions
source("./R_helpers/make_phyloseq_helper.R")
source("./R_helpers/analysis_helper.R")

# Select source file directory
project_directory <- "G:/Gedeelde drives/PRJ_MBAG/"
# Optionally convert Google Drive path to English
project_directory <- convert_gdrive_to_eng(project_directory, convert_to_eng = FALSE)

file_directory_bioinformatica <- paste0(project_directory, "4c_bodembiodiversiteit/data/bio-informatica")
file_directory_dataframe_overkoepelend <- paste0(project_directory, "4c_bodembiodiversiteit/data/statistiek/dataframe_overkoepelend/INBO+ILVO/Used\ phyloseq\ objects")
file_directory_data_ILVO_filtered <- paste0(project_directory, "4c_bodembiodiversiteit/data/ILVO")

# Path to Google Drive directory with source metadata file
source_data_dir <- "G:/Gedeelde drives/PRJ_MBAG/4c_bodembiodiversiteit/data"
# Path to Google Drive directory with source phyloseq objects
source_data_dir_bioinformatica <- file.path(source_data_dir, "bio-informatica")

# Load metadata
metadata_path <- paste0(project_directory, "4c_bodembiodiversiteit/data/Stratificatie_MBAG_plots/MBAG_stratfile_v2_cleaned_20.csv")
metadata <- load_metadata_mbag(metadata_path = metadata_path,
                               remove_mock_blancos_ntcs = TRUE)

# Set output directory
output_directory <- "./R_scripts/compare_primers"

# Define list of replicates that don't follow the standard replicate name system,
# and need to be separately addressed
custom_replicate_id_list <- c("137_012_241024_999_0579_207_01_25002",
                              "137_012_241024_999_0544_207_01_25002")

################################################################################

################################################################################

# 4. InseKP - All

# Make the full InseKP phyloseq object ("All") from the input files
physeq_InseKP_all <- make_insekp_24003_25002_phyloseq(only_samples=TRUE,
                                                      eng_gdrive=FALSE)

physeq_InseKP_all_bold <- physeq_InseKP_all

# Replace InseKP taxonomy file with boldigger3 db 3 mode 3 output

#see also script replace_taxonomy.R

# Load the new taxonomy file

taxonomy_df <- read.csv(file.path(source_data_dir_bioinformatica, "InseKP/data/24003_25002_fosf_combined/boldigger3_taxonomy/MBAG_InseKP_24003_25002_all_samples_MUMUFIED_boldigger3.csv"), stringsAsFactors = FALSE)

# Set Kingdom based on Phylum content
taxonomy_df$Kingdom <- ifelse(
  taxonomy_df$Phylum %in% c("unclassified_Root", "IncompleteTaxonomy"),
  "unclassified_Root",
  "Eukaryota"
)

# Define and reorder desired rank columns
original_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
lowercase_ranks <- tolower(original_ranks)

taxonomy_df <- taxonomy_df[, c("id", original_ranks)]

# Fill in unclassified levels
taxonomy_filled <- taxonomy_df
for (i in 1:nrow(taxonomy_df)) {
  last_classified <- NA
  for (col in original_ranks) {
    if (taxonomy_filled[i, col] == "" || is.na(taxonomy_filled[i, col])) {
      taxonomy_filled[i, col] <- paste0("unclassified_", last_classified)
    } else {
      last_classified <- taxonomy_filled[i, col]
    }
  }
}

# Rename columns to lowercase
colnames(taxonomy_filled)[2:8] <- lowercase_ranks  # Columns 2 to 8 are the taxonomy ranks

# Convert to matrix
taxonomy_mat <- as.matrix(taxonomy_filled[, lowercase_ranks])
rownames(taxonomy_mat) <- taxonomy_filled$id

# Replace taxonomy in phyloseq object
tax_table(physeq_InseKP_all_bold) <- tax_table(taxonomy_mat)

# ---

# InseKP - All

# Create phyloseq object with only known species
# Filter away entries with _otuXXX as suffix for unknown species
physeq_InseKP_all_bold_known_species <- subset_taxa(physeq_InseKP_all_bold, !grepl("_?unclassified_?|_otu\\d+", species))

# Agglomerate known species OTUs with tax_glom
physeq_InseKP_all_bold_known_species_tax_glom <- tax_glom(physeq_InseKP_all_bold_known_species, taxrank="species")

# Make diversity table for InseKP (full)
table_InseKP_full <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_all_bold,
  taxon_name = "insekp_full", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

# Make diversity table for InseKP (only known species)
table_InseKP_full_known_species_tax_glom <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_all_bold_known_species_tax_glom,
  taxon_name = "insekp_full_tax_glom", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

table_InseKP_full_known_species <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_all_bold_known_species,
  taxon_name = "insekp_full", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")


# InseKP - Annelida

# Subset the full InseKP to retain only phylum Annelida
physeq_InseKP_annelida = subset_taxa(physeq_InseKP_all, phylum=="Annelida")

# Create phyloseq object with only known species
# Filter away entries with _otuXXX as suffix for unknown species
physeq_InseKP_annelida_known_species <- subset_taxa(physeq_InseKP_annelida, !grepl("_?unclassified_?|_otu\\d+", species))

# Agglomerate known species OTUs with tax_glom
physeq_InseKP_annelida_known_species_tax_glom <- tax_glom(physeq_InseKP_annelida_known_species, taxrank="species")


# Make diversity table for InseKP Annelida (full)
table_InseKP_annelida <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_annelida,
  taxon_name = "annelida", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

# Make diversity table for InseKP Annelida (only known species)
table_InseKP_annelida_known_species_tax_glom <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_annelida_known_species_tax_glom,
  taxon_name = "annelida_tax_glom", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

table_InseKP_annelida_known_species <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_annelida_known_species,
  taxon_name = "annelida", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

# ---

# InseKP - Arthropoda

# Subset the full InseKP object to retain only phylum Arthropoda
physeq_InseKP_arthropoda = subset_taxa(physeq_InseKP_all, phylum=="Arthropoda")

# Create phyloseq object with only known species
# Filter away entries with _otuXXX as suffix for unknown species
physeq_InseKP_arthropoda_known_species <- subset_taxa(physeq_InseKP_arthropoda, !grepl("_?unclassified_?|_otu\\d+", species))

# Agglomerate known species OTUs with tax_glom
physeq_InseKP_arthropoda_known_species_tax_glom <- tax_glom(physeq_InseKP_arthropoda_known_species, taxrank="species")


# Make diversity table for InseKP Arthropoda (full)
table_InseKP_arthropoda <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_arthropoda,
  taxon_name = "arthropoda", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

# Make diversity table for InseKP Arthropoda (only known species)
table_InseKP_arthropoda_known_species_tax_glom <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_arthropoda_known_species_tax_glom,
  taxon_name = "arthropoda_tax_glom", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

table_InseKP_arthropoda_known_species <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_arthropoda_known_species,
  taxon_name = "arthropoda", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

# ---

# InseKP - Collembola

# Subset the full InseKP to retain only class Collembola
physeq_InseKP_collembola = subset_taxa(physeq_InseKP_all, class=="Collembola")

# Create phyloseq object with only known species
# Filter away entries with _otuXXX as suffix for unknown species
physeq_InseKP_collembola_known_species <- subset_taxa(physeq_InseKP_collembola, !grepl("_?unclassified_?|_otu\\d+", species))

# Agglomerate known species OTUs with tax_glom
physeq_InseKP_collembola_known_species_tax_glom <- tax_glom(physeq_InseKP_collembola_known_species, taxrank="species")


# Make diversity table for InseKP Collembola (full)
table_InseKP_collembola <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_collembola,
  taxon_name = "collembola", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

# Make diversity table for InseKP Collembola (only known species)
table_InseKP_collembola_known_species_tax_glom <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_collembola_known_species_tax_glom,
  taxon_name = "collembola_tax_glom", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

table_InseKP_collembola_known_species <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_collembola_known_species,
  taxon_name = "collembola", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

# ---

# Combined the InseKP tables vertically
combined_table_insekp <- bind_rows(
  table_InseKP_full,
  table_InseKP_annelida,
  table_InseKP_arthropoda,
  table_InseKP_collembola
)

# Print to file the InseKP combined tables
write_to_csv_and_rdata(diversity_index_table = combined_table_insekp,
                       output_directory = output_directory,
                       output_filename = "mbag_InseKP_dataframe.csv")


# Combined the InseKP known_species_tax_glom tables vertically
combined_table_insekp_known_species_tax_glom <- bind_rows(
  table_InseKP_full_known_species_tax_glom,
  table_InseKP_annelida_known_species_tax_glom,
  table_InseKP_arthropoda_known_species_tax_glom,
  table_InseKP_collembola_known_species_tax_glom,
  table_InseKP_full_known_species,
  table_InseKP_annelida_known_species,
  table_InseKP_arthropoda_known_species,
  table_InseKP_collembola_known_species
)

# Print to file the InseKP known_species_tax_glom combined tables
write_to_csv_and_rdata(diversity_index_table = combined_table_insekp_known_species_tax_glom,
                       output_directory = output_directory,
                       output_filename = "mbag_InseKP_dataframe_known_species.csv")

#### InseKP based on BOLD data

# ---
# InseKP - Annelida (BOLD taxonomy)

# Subset the InseKP object with BOLD taxonomy to retain only phylum Annelida
physeq_InseKP_annelida_bold = subset_taxa(physeq_InseKP_all_bold, phylum=="Annelida")

# Filter away entries with _otuXXX as suffix or unclassified species
physeq_InseKP_annelida_bold_known_species <- subset_taxa(physeq_InseKP_annelida_bold, !grepl("_?unclassified_?|_otu\\d+", species))

# Agglomerate known species OTUs with tax_glom
physeq_InseKP_annelida_bold_known_species_tax_glom <- tax_glom(physeq_InseKP_annelida_bold_known_species, taxrank="species")

# Make diversity table (full)
table_InseKP_annelida_bold <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_annelida_bold,
  taxon_name = "annelida_bold", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

# Make diversity table (only known species)
table_InseKP_annelida_bold_known_species_tax_glom <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_annelida_bold_known_species_tax_glom,
  taxon_name = "annelida_bold_tax_glom", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

table_InseKP_annelida_bold_known_species <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_annelida_bold_known_species,
  taxon_name = "annelida_bold", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

# ---
# InseKP - Lumbricidae (BOLD taxonomy)

# Subset the InseKP object with BOLD taxonomy to retain only phylum Annelida
physeq_InseKP_Lumbricidae_bold = subset_taxa(physeq_InseKP_annelida_bold, family=="Lumbricidae")

# Filter away entries with _otuXXX as suffix or unclassified species
physeq_InseKP_Lumbricidae_bold_known_species <- subset_taxa(physeq_InseKP_Lumbricidae_bold, !grepl("_?unclassified_?|_otu\\d+", species))

# Agglomerate known species OTUs with tax_glom
physeq_InseKP_Lumbricidae_bold_known_species_tax_glom <- tax_glom(physeq_InseKP_Lumbricidae_bold_known_species, taxrank="species")

# Make diversity table (full)
table_InseKP_lumbricidae_bold <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_Lumbricidae_bold,
  taxon_name = "lumbricidae_bold", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

# Make diversity table (only known species)
table_InseKP_lumbricidae_bold_known_species_tax_glom <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_Lumbricidae_bold_known_species_tax_glom,
  taxon_name = "lumbricidae_bold_tax_glom", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

table_InseKP_lumbricidae_bold_known_species <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_Lumbricidae_bold_known_species,
  taxon_name = "lumbricidae_bold", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

# ---
# InseKP - Enchytraeidae (BOLD taxonomy)

# Subset the InseKP object with BOLD taxonomy to retain only phylum Annelida
physeq_InseKP_Enchytraeidae_bold = subset_taxa(physeq_InseKP_annelida_bold, family=="Enchytraeidae")

# Filter away entries with _otuXXX as suffix or unclassified species
physeq_InseKP_Enchytraeidae_bold_known_species <- subset_taxa(physeq_InseKP_Enchytraeidae_bold, !grepl("_?unclassified_?|_otu\\d+", species))

# Agglomerate known species OTUs with tax_glom
physeq_InseKP_Enchytraeidae_bold_known_species_tax_glom <- tax_glom(physeq_InseKP_Enchytraeidae_bold_known_species, taxrank="species")

# Make diversity table (full)
table_InseKP_enchytraeidae_bold <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_Enchytraeidae_bold,
  taxon_name = "enchytraeidae_bold", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

# Make diversity table (only known species)
table_InseKP_enchytraeidae_bold_known_species_tax_glom <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_Enchytraeidae_bold_known_species_tax_glom,
  taxon_name = "enchytraeidae_bold_tax_glom", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

table_InseKP_enchytraeidae_bold_known_species <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_Enchytraeidae_bold_known_species,
  taxon_name = "enchytraeidae_bold", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

# ---
# InseKP - Arthropoda (BOLD taxonomy)

physeq_InseKP_arthropoda_bold = subset_taxa(physeq_InseKP_all_bold, phylum=="Arthropoda")
physeq_InseKP_arthropoda_bold_known_species <- subset_taxa(physeq_InseKP_arthropoda_bold, !grepl("_?unclassified_?|_otu\\d+", species))
physeq_InseKP_arthropoda_bold_known_species_tax_glom <- tax_glom(physeq_InseKP_arthropoda_bold_known_species, taxrank="species")

table_InseKP_arthropoda_bold <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_arthropoda_bold,
  taxon_name = "arthropoda_bold", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

table_InseKP_arthropoda_bold_known_species_tax_glom <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_arthropoda_bold_known_species_tax_glom,
  taxon_name = "arthropoda_bold_tax_glom", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

table_InseKP_arthropoda_bold_known_species <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_arthropoda_bold_known_species,
  taxon_name = "arthropoda_bold", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

# ---
# InseKP - Collembola (BOLD taxonomy)

physeq_InseKP_collembola_bold = subset_taxa(physeq_InseKP_all_bold, class=="Collembola")
physeq_InseKP_collembola_bold_known_species <- subset_taxa(physeq_InseKP_collembola_bold, !grepl("_?unclassified_?|_otu\\d+", species))
physeq_InseKP_collembola_bold_known_species_tax_glom <- tax_glom(physeq_InseKP_collembola_bold_known_species, taxrank="species")

table_InseKP_collembola_bold <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_collembola_bold,
  taxon_name = "collembola_bold", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

table_InseKP_collembola_bold_known_species_tax_glom <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_collembola_bold_known_species_tax_glom,
  taxon_name = "collembola_bold_tax_glom", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

table_InseKP_collembola_bold_known_species <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_InseKP_collembola_bold_known_species,
  taxon_name = "collembola_bold", unit_value="swarm", primerset = "InseKP",
  metadata_df = metadata,
  custom_replicate_id_list = custom_replicate_id_list,
  micro_or_macro_organism = "MACRO")

# ---

# Combine BOLD tables with the original ones
combined_table_insekp <- bind_rows(
  combined_table_insekp,
  table_InseKP_annelida_bold,
  table_InseKP_arthropoda_bold,
  table_InseKP_collembola_bold,
  table_InseKP_enchytraeidae_bold,
  table_InseKP_lumbricidae_bold
)

combined_table_insekp_known_species_tax_glom <- bind_rows(
  combined_table_insekp_known_species_tax_glom,
  table_InseKP_annelida_bold_known_species_tax_glom,
  table_InseKP_arthropoda_bold_known_species_tax_glom,
  table_InseKP_collembola_bold_known_species_tax_glom,
  table_InseKP_enchytraeidae_bold_known_species_tax_glom,
  table_InseKP_lumbricidae_bold_known_species_tax_glom,
  table_InseKP_annelida_bold_known_species,
  table_InseKP_lumbricidae_bold_known_species,
  table_InseKP_enchytraeidae_bold_known_species,
  table_InseKP_arthropoda_bold_known_species,
  table_InseKP_collembola_bold_known_species
)

# Save updated combined tables including BOLD
write_to_csv_and_rdata(diversity_index_table = combined_table_insekp,
                       output_directory = output_directory,
                       output_filename = "mbag_InseKP_dataframe_incl_BOLD.csv")

write_to_csv_and_rdata(diversity_index_table = combined_table_insekp_known_species_tax_glom,
                       output_directory = output_directory,
                       output_filename = "mbag_InseKP_dataframe_known_species_incl_BOLD.csv")



####
# Add total_total_count (all)

# Get total counts per sample from the full phyloseq object
total_counts_Insekp_all <- data.frame(
  sample = sample_names(physeq_InseKP_all),
  total_total_count = sample_sums(physeq_InseKP_all)
)

# All (with unknown species)
# Merge combined diversity index table with total counts
combined_table_InseKP_INBO_with_total <- combined_table_insekp %>%
  dplyr::left_join(total_counts_Insekp_all, by = "sample")

# Save the updated combined table
write_to_csv_and_rdata(diversity_index_table = combined_table_InseKP_INBO_with_total,
                       output_directory = output_directory,
                       output_filename = "mbag_InseKP_INBO_dataframe_with_total_counts.csv")

# All (only known species)
# Merge combined diversity index table with total counts
combined_table_InseKP_INBO_with_total_known_species <- combined_table_insekp_known_species_tax_glom %>%
  dplyr::left_join(total_counts_Insekp_all, by = "sample")

# Save the updated combined table
write_to_csv_and_rdata(diversity_index_table = combined_table_InseKP_INBO_with_total_known_species,
                       output_directory = output_directory,
                       output_filename = "mbag_InseKP_INBO_dataframe_with_total_counts_known_species.csv")



################################################################################

# Include BOLD-InseKP in full (with unknown species) INBO combined table
combined_table_mbag_INBO <- bind_rows(combined_table_18S_INBO,
                                      table_Coll01_collembola,
                                      table_Olig01_annelida,
                                      combined_table_insekp)

# Write updated INBO combined table
write_to_csv_and_rdata(diversity_index_table = combined_table_mbag_INBO,
                       output_directory = output_directory,
                       output_filename = "mbag_combined_dataframe_INBO.csv")

# Include BOLD-InseKP in known species INBO combined table
combined_table_mbag_INBO_known_species <- bind_rows(table_Coll01_collembola_known_species_tax_glom,
                                                    table_Coll01_collembola_known_species,
                                                    table_Olig01_annelida_known_species_tax_glom,
                                                    table_Olig01_annelida_known_species,
                                                    combined_table_insekp_known_species_tax_glom)

# Write updated INBO known-species combined table
write_to_csv_and_rdata(diversity_index_table = combined_table_mbag_INBO_known_species,
                       output_directory = output_directory,
                       output_filename = "mbag_combined_dataframe_INBO_known_species.csv")


####
# Include BOLD-InseKP total-count table in full (with unknown) INBO combined table
combined_table_mbag_INBO_all <- bind_rows(combined_table_18S_INBO_with_total,
                                          combined_table_Coll01_INBO_with_total,
                                          combined_table_Olig01_INBO_with_total,
                                          combined_table_InseKP_INBO_with_total)

write_to_csv_and_rdata(diversity_index_table = combined_table_mbag_INBO_all,
                       output_directory = output_directory,
                       output_filename = "mbag_combined_dataframe_INBO_all.csv")

# Include BOLD-InseKP total-count table in known-species INBO combined table
combined_table_mbag_INBO_known_species_all <- bind_rows(combined_table_Coll01_INBO_with_total_known_species,
                                                        combined_table_Coll01_INBO_with_total_known_species_tax_glom,
                                                        combined_table_Olig01_INBO_with_total_known_species,
                                                        combined_table_Olig01_INBO_with_total_known_species_tax_glom,
                                                        combined_table_InseKP_INBO_with_total_known_species)

write_to_csv_and_rdata(diversity_index_table = combined_table_mbag_INBO_known_species_all,
                       output_directory = output_directory,
                       output_filename = "mbag_combined_dataframe_INBO_known_species_total_counts.csv")





#-----------
# Add total_count (full dataset)

# Combined the Fungi tables vertically
combined_table_fungi_ilvo <- bind_rows(
  table_ILVO_Fungi,
  table_ILVO_Fungi_filtered
)

# Print to file the Fungi combined tables
write_to_csv_and_rdata(diversity_index_table = combined_table_fungi_ilvo,
                       output_directory = output_directory,
                       output_filename = "mbag_fungi_dataframe_ilvo.csv")


# Add total_total_count (all)

# Get total counts per sample from the full phyloseq object
total_counts_fungi_ilvo_all <- data.frame(
  sample = sample_data(physeq_ILVO_Fungi)$Sample,
  total_total_count = sample_sums(physeq_ILVO_Fungi)
)

# All (with unknown species)
# Merge combined diversity index table with total counts
combined_table_fungi_ILVO_with_total <- combined_table_fungi_ilvo %>%
  dplyr::left_join(total_counts_fungi_ilvo_all, by = "sample")

# Save the updated combined table
write_to_csv_and_rdata(diversity_index_table = combined_table_fungi_ILVO_with_total,
                       output_directory = output_directory,
                       output_filename = "mbag_fungi_ILVO_dataframe_with_total_counts.csv")





# All (with unknown species)
# Merge combined diversity index table with total counts
combined_table_fungi_ILVO_phyla_with_total <- combined_all_diversity_tables_phyla_fungi %>%
  dplyr::left_join(total_counts_fungi_ilvo_phyla_all, by = "sample")


# observed is currently NA, should have observed = 0, because the data in this case is compositional
combined_table_fungi_ILVO_phyla_with_total <- combined_table_fungi_ILVO_phyla_with_total %>%
  mutate(
    observed = as.numeric(observed),
    observed = if_else(is.na(observed), 0, observed)
  )


#-----------

# 5. Nematoden

# Load in the non-rarefied dataframe
physeq_ILVO_Nem_NA <- readRDS(file.path(file_directory_dataframe_overkoepelend, "MBAG_18S_with_metadata.rds"))
# Bring sample names to lower-case
sample_data(physeq_ILVO_Nem_NA)$Sample <- tolower(sample_data(physeq_ILVO_Nem_NA)$Sample)

# ---

# Make physeq_ILVO_Nem_filtered
#
# Subset to remove ASVs with undefined genus level
# Note: The taxonomy annotation contains "-NA" at genus level
#       (e.g. "ASVXXX-NA") and their row name is the same.

# Extract the tax_table and otu_table
tax_table_df <- as.data.frame(phyloseq::tax_table(physeq_ILVO_Nem_NA))
otu_table_df <- as.data.frame(phyloseq::otu_table(physeq_ILVO_Nem_NA))

# Identify row names containing "-NA" in tax_table and otu_table
# Note: grep is by default case-sensitive, so the search string "-NA" is safe
rows_to_remove <- grep("-NA", rownames(tax_table_df))

# Remove rows from tax_table
if (length(rows_to_remove) > 0) {
  tax_table_filtered <- tax_table_df[-rows_to_remove, , drop = FALSE]
} else {
  tax_table_filtered <- tax_table_df # No rows to remove
}

# Remove rows from otu_table
rows_to_remove_otu <- grep("-NA", rownames(otu_table_df))

if (length(rows_to_remove_otu) > 0) {
  otu_table_filtered <- otu_table_df[-rows_to_remove_otu, , drop = FALSE]
} else {
  otu_table_filtered <- otu_table_df # No rows to remove
}

# Update the phyloseq object with the filtered tax_table and otu_table
physeq_ILVO_Nem <- phyloseq::phyloseq(
  phyloseq::otu_table(otu_table_filtered, taxa_are_rows = TRUE),
  phyloseq::tax_table(as.matrix(tax_table_filtered)),
  phyloseq::sample_data(phyloseq::sample_data(physeq_ILVO_Nem_NA))
)

# ---

# Subset to data with Material "eDNA" and "nema-extract"
physeq_ILVO_Nem_eDNA <- subset_samples(physeq_ILVO_Nem, Material == "eDNA")
physeq_ILVO_Nem_nema_extract <- subset_samples(physeq_ILVO_Nem, Material == "nema-extract")

# ---

# Make diversity table for Nematodes "eDNA"
table_ILVO_Nem_eDNA <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_ILVO_Nem_eDNA,
  taxon_name = "Nematoda", unit_value = "ASV", primerset = "18S_nema_edna",
  metadata_df = metadata,
  micro_or_macro_organism = "micro",
  tidytacos_sample_header = "Sample",
  list_runs_with_unrelated_samples = c("25002"),
  add_missing_samples = FALSE)


# Make diversity table for Nematodes "nema extract"
table_ILVO_Nem_nema_extract  <- make_diversity_index_table_from_phyloseq_obj(
  phyloseq_obj = physeq_ILVO_Nem_nema_extract,
  taxon_name = "Nematoda", unit_value = "ASV", primerset = "18S_nema_extract",
  metadata_df = metadata,
  micro_or_macro_organism = "micro",
  tidytacos_sample_header = "Sample",
  list_runs_with_unrelated_samples = c("25002"),
  add_missing_samples = FALSE)


#-----------
# Add total_count (full dataset)

# Combined the InseKP tables vertically
combined_table_18S_nem_ilvo <- bind_rows(
  table_ILVO_Nem_eDNA,
  table_ILVO_Nem_nema_extract
)

# Print to file the InseKP combined tables
write_to_csv_and_rdata(diversity_index_table = combined_table_18S_nem_ilvo,
                       output_directory = output_directory,
                       output_filename = "mbag_18S_nem_dataframe_ilvo.csv")


# Add total_total_count (all)

# Get total counts per sample from the full phyloseq object
total_counts_18S_nem_ilvo_all <- data.frame(
  sample = sample_data(physeq_ILVO_Nem)$Sample,
  total_total_count = sample_sums(physeq_ILVO_Nem)
)

# All (with unknown species)
# Merge combined diversity index table with total counts
combined_table_18S_nem_ILVO_with_total <- combined_table_18S_nem_ilvo %>%
  dplyr::left_join(total_counts_fungi_ilvo_all, by = "sample")

# Save the updated combined table
write_to_csv_and_rdata(diversity_index_table = combined_table_18S_nem_ILVO_with_total,
                       output_directory = output_directory,
                       output_filename = "mbag_18S_nem_ILVO_dataframe_with_total_counts.csv")

#--------
#--------

# Combine the diversity index tables for all ILVO primer sets into one table

# Helper function to convert 'observed' to numeric
ensure_numeric_observed <- function(df) {
  df %>% mutate(observed = as.numeric(observed))
}

# Apply to all data frames
combined_table_mbag_ILVO_all <- bind_rows(
  ensure_numeric_observed(combined_table_16S_ILVO_with_total),
  ensure_numeric_observed(combined_table_16S_ILVO_phyla_with_total),
  ensure_numeric_observed(combined_table_fungi_ILVO_with_total),
  ensure_numeric_observed(combined_table_fungi_ILVO_phyla_with_total),
  ensure_numeric_observed(combined_table_18S_nem_ILVO_with_total)
)


write_to_csv_and_rdata(diversity_index_table = combined_table_mbag_ILVO_all,
                       output_directory = output_directory,
                       output_filename = "mbag_ILVO_dataframe.csv")


#--------
#--------
# Total_counts all full dataset

# Combine the diversity index tables for all ILVO primer sets into one table

combined_table_mbag_ILVO_all <- bind_rows(
  ensure_numeric_observed(combined_table_16S_ILVO_with_total),
  ensure_numeric_observed(combined_table_16S_ILVO_phyla_with_total),
  ensure_numeric_observed(combined_table_fungi_ILVO_with_total),
  ensure_numeric_observed(combined_table_fungi_ILVO_phyla_with_total),
  ensure_numeric_observed(combined_table_18S_nem_ILVO_with_total)
)

write_to_csv_and_rdata(diversity_index_table = combined_table_mbag_ILVO_all,
                       output_directory = output_directory,
                       output_filename = "mbag_ILVO_dataframe_all_total_counts.csv")

###############################################################################

# Combine the diversity index tables for all INBO and ILVO primer sets into one table

combined_table_mbag_INBO$observed <- as.numeric(combined_table_mbag_INBO$observed)
combined_table_mbag_ILVO_all$observed <- as.numeric(combined_table_mbag_ILVO_all$observed)


# Stack all tables vertically
combined_table_INBO_ILVO <- bind_rows(
  combined_table_mbag_INBO,
  combined_table_mbag_ILVO_all
)


# Unrelated to the bacterial and fungal phyla, but in general for al entries in the genetic and taxonomic combined dataframes, when observed is 0, NA values for shannon and simpson should be replaced by 0

# Ensure columns are numeric, then replace NA with 0 if observed == 0
combined_table_INBO_ILVO <- combined_table_INBO_ILVO %>%
  mutate(
    shannon = as.numeric(shannon),
    simpson = as.numeric(simpson),
    shannon = if_else(observed == 0 & is.na(shannon), 0, shannon),
    simpson = if_else(observed == 0 & is.na(simpson), 0, simpson)
  )




# Write csv
write_to_csv_and_rdata(diversity_index_table = combined_table_INBO_ILVO,
                       output_directory = output_directory,
                       output_filename = "mbag_combined_dataframe_INBO_ILVO.csv")


##--------
# All total_counts

# Stack all tables vertically
combined_table_INBO_ILVO_total_counts <- bind_rows(
  ensure_numeric_observed(combined_table_mbag_INBO_all),
  ensure_numeric_observed(combined_table_mbag_ILVO_all)
)

# Unrelated to the bacterial and fungal phyla, but in general for al entries in the genetic and taxonomic combined dataframes, when observed is 0, NA values for shannon and simpson should be replaced by 0

# Ensure columns are numeric, then replace NA with 0 if observed == 0
combined_table_INBO_ILVO_total_counts <- combined_table_INBO_ILVO_total_counts %>%
  mutate(
    shannon = as.numeric(shannon),
    simpson = as.numeric(simpson),
    shannon = if_else(observed == 0 & is.na(shannon), 0, shannon),
    simpson = if_else(observed == 0 & is.na(simpson), 0, simpson)
  )

write_to_csv_and_rdata(diversity_index_table = combined_table_INBO_ILVO_total_counts,
                       output_directory = output_directory,
                       output_filename = "mbag_combined_dataframe_INBO_ILVO_total_counts.csv")

# Remove samples for which the eDNA metabarcoding protocol failed

# Filter out rows where 'primerset' is "InseKP" and 'sample' matches the specified names
combined_table_INBO_ILVO_total_counts <- combined_table_INBO_ILVO_total_counts %>%
  filter(!(primerset == "InseKP" & sample %in% c("137_012_241022_999_0438_207_01_25002", "P89ecc44_0_10")))


# Write csv
write_to_csv_and_rdata(diversity_index_table = combined_table_INBO_ILVO_total_counts,
                       output_directory = output_directory,
                       output_filename = "mbag_combined_dataframe_INBO_ILVO_total_counts_final.csv")
################################################################################

# Add metadata to each INBO+ILVO sample rows in the combined dataframes

# All (including unknown species)

combined_table_INBO_ILVO_metadata <- merge_combined_diversity_table_with_metadata(
  combined_diversity_index_table = combined_table_INBO_ILVO,
  metadata_df = metadata,
  combined_table_sample_column_header = "sample",
  metadata_df_sample_column_header = "row.names")

# Save
write_to_csv_and_rdata(diversity_index_table = combined_table_INBO_ILVO_metadata,
                       output_directory = output_directory,
                       output_filename = "mbag_combined_dataframe_INBO_ILVO_metadata.csv")

# ---

# Only known species

combined_table_mbag_INBO_known_species_metadata <- merge_combined_diversity_table_with_metadata(
  combined_diversity_index_table = combined_table_mbag_INBO_known_species,
  metadata_df = metadata,
  combined_table_sample_column_header = "sample",
  metadata_df_sample_column_header = "row.names")
# Save
write_to_csv_and_rdata(diversity_index_table = combined_table_mbag_INBO_known_species_metadata,
                       output_directory = output_directory,
                       output_filename = "mbag_combined_dataframe_INBO_known_species_metadata.csv")



##### -----------
# All total_counts + metadata

combined_table_INBO_ILVO_total_counts_metadata <- merge_combined_diversity_table_with_metadata(
  combined_diversity_index_table = combined_table_INBO_ILVO_total_counts,
  metadata_df = metadata,
  combined_table_sample_column_header = "sample",
  metadata_df_sample_column_header = "row.names")

# Save
write_to_csv_and_rdata(diversity_index_table = combined_table_INBO_ILVO_metadata,
                       output_directory = output_directory,
                       output_filename = "mbag_combined_dataframe_INBO_ILVO_total_counts_metadata.csv")



