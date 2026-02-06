# Patch the R script to normalize OTU/feature IDs so that taxonomy 'swarm' like 1,2,3
# will match OTU IDs like 'SWARM1','SWARM2','SWARM3'. We standardize BOTH sides to 'SWARM<digits>'.
# Robust loader to build a phyloseq object from:
# - OTU table (tab-delimited): dieetanalyse_weidevogels_otu_table.table
# - Sample metadata (CSV):    dieetanalyse_weidevogels_metadata_curated.csv
# - Taxonomy table (CSV):     dieetanalyse_weidevogels_taxonomy_curated.csv
#
# This script attempts to auto-detect ID columns and clean formats to be compatible with phyloseq.

suppressPackageStartupMessages({
  library(phyloseq)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(magrittr)
})

# ---- File paths (edit if needed) ----
otu_path <- "25004_InseKP_SeDNA_EJPsoil_otu_table.table"
meta_path <- "metadata_testpooling.csv"
tax_path <- "25004_InseKP_SeDNA_EJPsoil_crabs_e02_tax.csv"

# point to the current working directory
otu_path  <- file.path(getwd(), basename(otu_path))
meta_path <- file.path(getwd(), basename(meta_path))
tax_path  <- file.path(getwd(), basename(tax_path))

message("Reading files:\n - OTU: ", otu_path,
        "\n - META: ", meta_path,
        "\n - TAX: ", tax_path)


# Helper: normalize any feature ID to canonical 'SWARM<digits>' if it looks like a swarm id
norm_swarm_id <- function(x) {
  x <- as.character(x)
  # if already like SWARM123 (any case), standardize case
  already <- str_match(x, "^[sS][wW][aA][rR][mM]\\s*0*([0-9]+)$")[,2]
  onlydigits <- str_match(x, "^\\s*0*([0-9]+)\\s*$")[,2]
  out <- ifelse(!is.na(already),
                paste0("SWARM", already),
                ifelse(!is.na(onlydigits),
                       paste0("SWARM", onlydigits),
                       x))
  out
}

# ---- Read OTU table ----
otu_raw <- read_delim(otu_path, delim = "\t", col_types = cols(.default = col_guess()))
otu_id_col <- names(otu_raw)[1]
cand_id_cols <- c("#otu id","swarm","otu","otu_id","featureid","feature","asv","asv_id","id")
id_hit <- intersect(tolower(names(otu_raw)), cand_id_cols)
if (length(id_hit) == 1) {
  otu_id_col <- names(otu_raw)[match(id_hit, tolower(names(otu_raw)))]
}
stopifnot(otu_id_col %in% names(otu_raw))

otu_tab <- otu_raw %>%
  rename(.feature_id = all_of(otu_id_col)) %>%
  mutate(.feature_id = norm_swarm_id(.feature_id)) %>%
  distinct(.feature_id, .keep_all = TRUE) %>%
  column_to_rownames(".feature_id")
otu_tab[] <- lapply(otu_tab, function(x) suppressWarnings(as.numeric(x)))
otu_tab[is.na(otu_tab)] <- 0
taxa_are_rows_flag <- TRUE

# ---- Read sample metadata ----
meta_raw <- read_csv(meta_path, show_col_types = FALSE)
cand_sid <- c("sample","sample_id","sampleid","samplename","library","run","sample_name","id")
sid_hit <- intersect(tolower(names(meta_raw)), cand_sid)
if (length(sid_hit) >= 1) {
  sample_id_col <- names(meta_raw)[match(sid_hit[1], tolower(names(meta_raw)))]
} else {
  sample_id_col <- names(meta_raw)[1]
  message("No obvious sample ID column found in metadata; using first column: ", sample_id_col)
}
meta_tbl <- meta_raw %>%
  rename(.sample_id = all_of(sample_id_col)) %>%
  mutate(.sample_id = as.character(.sample_id)) %>%
  distinct(.sample_id, .keep_all = TRUE) %>%
  column_to_rownames(".sample_id")
meta_tbl <- meta_tbl %>%
  mutate(across(where(is.logical), as.character)) %>%
  mutate(across(where(is.list), ~vapply(., function(x) paste(x, collapse=";"), character(1))))

# ---- Read taxonomy table ----
tax_raw <- read_csv(tax_path, show_col_types = FALSE)
tax_id_col <- names(tax_raw)[1]
id_hit2 <- intersect(tolower(names(tax_raw)), c("swarm", "swarm_id", "otu","otu_id","featureid","feature","asv","asv_id","id"))
if (length(id_hit2) == 1) {
  tax_id_col <- names(tax_raw)[match(id_hit2, tolower(names(tax_raw)))]
}
tax_tbl <- tax_raw %>%
  rename(.feature_id = all_of(tax_id_col)) %>%
  mutate(.feature_id = norm_swarm_id(.feature_id)) %>%
  distinct(.feature_id, .keep_all = TRUE)

rank_map <- c(
  "kingdom"="kingdom",
  "phylum"="phylum",
  "class"="class",
  "order"="order",
  "family"="family",
  "genus"="genus",
  "species"="species"
)
tax_cols <- lapply(names(rank_map), function(rk) {
  ix <- match(rk, tolower(names(tax_tbl)))
  if (!is.na(ix)) names(tax_tbl)[ix] else NA_character_
}) %>% setNames(names(rank_map))
tax_rank_cols <- unname(na.omit(unlist(tax_cols)))
if (length(tax_rank_cols) == 0) {
  stop("No recognizable taxonomic rank columns (kingdom..species) found in taxonomy file.")
}
tax_tbl <- tax_tbl %>%
  mutate(across(all_of(tax_rank_cols), ~na_if(str_trim(as.character(.x)), "")))

# Capitalize selected ranks: phylum, class, order, family, genus
capitalize_ranks <- c("phylum","class","order","family","genus")
for (rk in capitalize_ranks) {
  if (rk %in% tolower(names(tax_tbl))) {
    rk_col <- names(tax_tbl)[match(rk, tolower(names(tax_tbl)))]
    tax_tbl[[rk_col]] <- ifelse(is.na(tax_tbl[[rk_col]]), NA_character_,
                                str_replace_all(str_to_sentence(tax_tbl[[rk_col]]), "_", " "))
  }
}

tax_mat <- tax_tbl %>%
  select(.feature_id, all_of(tax_rank_cols)) %>%
  column_to_rownames(".feature_id") %>%
  as.matrix()

# ---- Align keys across tables ----
shared_samples <- intersect(colnames(otu_tab), rownames(meta_tbl))
if (length(shared_samples) == 0) {
  stop("No overlapping sample IDs between OTU table (columns) and metadata (rownames).")
}
if (length(shared_samples) < ncol(otu_tab)) {
  message("Dropping ", ncol(otu_tab) - length(shared_samples),
          " samples from OTU table not present in metadata.")
}
if (length(shared_samples) < nrow(meta_tbl)) {
  message("Dropping ", nrow(meta_tbl) - length(shared_samples),
          " metadata rows not present in OTU table.")
}
otu_tab <- otu_tab[, shared_samples, drop = FALSE]
meta_tbl <- meta_tbl[shared_samples, , drop = FALSE]

shared_features <- intersect(rownames(otu_tab), rownames(tax_mat))
if (length(shared_features) == 0) {
  # Provide diagnostics to help the user
  ots <- head(rownames(otu_tab), 10)
  txs <- head(rownames(tax_mat), 10)
  stop(paste0("No overlapping feature/OTU IDs between OTU table (rows) and taxonomy (rownames).\n",
              "Example OTU IDs: ", paste(ots, collapse=", "), "\n",
              "Example TAX IDs: ", paste(txs, collapse=", "), "\n",
              "Hint: IDs are normalized to 'SWARM<digits>'. Check that both sides follow this pattern."))
}
if (length(shared_features) < nrow(otu_tab)) {
  message("Dropping ", nrow(otu_tab) - length(shared_features),
          " OTUs from OTU table without taxonomy.")
}
if (length(shared_features) < nrow(tax_mat)) {
  message("Dropping ", nrow(tax_mat) - length(shared_features),
          " taxonomy rows without counts.")
}
otu_tab <- otu_tab[shared_features, , drop = FALSE]
tax_mat <- tax_mat[shared_features, , drop = FALSE]

# ---- Construct phyloseq object ----
ps <- phyloseq(
  otu_table(as.matrix(otu_tab), taxa_are_rows = TRUE),
  sample_data(as.data.frame(meta_tbl)),
  tax_table(tax_mat)
)

# ---- Report ----
message("phyloseq object created.")
message("  Samples: ", nsamples(ps))
message("  Taxa:    ", ntaxa(ps))
message("  Ranks:   ", paste(rank_names(ps), collapse = ", "))

print(ps)

saveRDS(ps, file = "phyloseq_object.rds")

