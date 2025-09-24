# sample IDs present in OTU table but missing from metadata
samples_in_otu_not_meta <- setdiff(colnames(otu_tab), rownames(meta_tbl))
length(samples_in_otu_not_meta); head(samples_in_otu_not_meta)

# sample IDs present in metadata but missing from OTU table
samples_in_meta_not_otu <- setdiff(rownames(meta_tbl), colnames(otu_tab))
length(samples_in_meta_not_otu); head(samples_in_meta_not_otu)

# inspect the actual entries that would be dropped
otu_missing_samples  <- otu_tab[, colnames(otu_tab) %in% samples_in_otu_not_meta, drop = FALSE]
meta_missing_samples <- meta_tbl[rownames(meta_tbl) %in% samples_in_meta_not_otu, , drop = FALSE]
