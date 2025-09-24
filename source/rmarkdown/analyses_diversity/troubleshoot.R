# IDs present in OTU table but NOT in metadata (these columns get dropped)
otu_only <- setdiff(colnames(otu_tab), shared_samples)
length(otu_only); head(otu_only)

# IDs present in metadata but NOT in OTU table (these rows get dropped)
meta_only <- setdiff(rownames(meta_tbl), shared_samples)
length(meta_only); head(meta_only)

# Inspect the actual entries that will be dropped
otu_missing  <- otu_tab[, otu_only, drop = FALSE]        # columns from OTU table
meta_missing <- meta_tbl[meta_only, , drop = FALSE]      # rows from metadata
