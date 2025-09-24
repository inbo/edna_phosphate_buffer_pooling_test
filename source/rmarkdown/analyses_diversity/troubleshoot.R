# --- Side-by-side list of sample IDs ---
otu_ids  <- colnames(otu_tab)
meta_ids <- rownames(meta_tbl)

# pad vectors to equal length with NA
pad <- function(x, n) { length(x) <- n; x }
n_max <- max(length(otu_ids), length(meta_ids))

side_by_side <- data.frame(
  otu_tab_colname   = pad(otu_ids,  n_max),
  meta_tbl_rownames = pad(meta_ids, n_max),
  stringsAsFactors = FALSE
)

write.csv(side_by_side, "sample_ids_side_by_side.csv", row.names = FALSE)


all_ids <- sort(unique(c(otu_ids, meta_ids)))
recon <- data.frame(
  sample_id = all_ids,
  in_otu  = all_ids %in% otu_ids,
  in_meta = all_ids %in% meta_ids,
  stringsAsFactors = FALSE
)
recon$matched_both <- recon$in_otu & recon$in_meta

write.csv(recon, "sample_id_reconciliation.csv", row.names = FALSE)

