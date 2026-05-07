# Spike-in identification and adjustment for phyloseq objects.
#
# Identifies spike-in ASVs using edit-distance (adist) fuzzy matching against
# in-silico V3V4 amplicons extracted from full-length 16S references. This
# replaces the previous tidysq::%has% exact substring approach, which missed
# ASVs with small sequence differences caused by binned quality scores and
# DADA2's primer/spacer trimming.
#
# Adds spike-in relative abundance and counts to sample_data, then removes
# spike-in ASVs from the phyloseq OTU table.

library(phyloseq)
library(microViz)
library(Biostrings)
library(tidyverse)

ps <- readRDS(snakemake@input[["ps"]])

# Load spike-in full-length 16S references
allobacillus_fasta <- readDNAStringSet(snakemake@input[["allobacillus_fasta"]])
imtechella_fasta <- readDNAStringSet(snakemake@input[["imtechella_fasta"]])
manual_taxa <- read_csv(snakemake@input[["spikein_manual_taxonomy"]])
spiked_in_species <- c("spike_in_Allobacillus_halotolerans", "spike_in_Imtechella_halotolerans")

# Primers and trimming params from Snakemake rule
fwd_primer <- snakemake@params[["forward_primer"]]
rev_primer <- snakemake@params[["reverse_primer"]]
trim_left <- as.integer(snakemake@params[["trim_left"]])
trim_right <- as.integer(snakemake@params[["trim_right"]])
max_edit_dist <- as.integer(snakemake@params[["max_edit_dist"]])

# --- Extract V3V4 amplicon from full-length 16S via in-silico PCR ---
# Finds primer binding sites, extracts the region between them, then trims
# the same number of bases that DADA2's trimLeft removes from each end.
extract_v3v4 <- function(full_16s, fwd_primer, rev_primer,
                         trim_left = 0, trim_right = 0, max_mm = 2) {
  fwd_pattern <- DNAString(fwd_primer)
  rev_pattern <- reverseComplement(DNAString(rev_primer))

  results <- character(0)
  for (i in seq_along(full_16s)) {
    seq <- full_16s[[i]]
    fwd_hits <- matchPattern(fwd_pattern, seq, max.mismatch = max_mm, fixed = "subject")
    rev_hits <- matchPattern(rev_pattern, seq, max.mismatch = max_mm, fixed = "subject")

    if (length(fwd_hits) > 0 && length(rev_hits) > 0) {
      amp_start <- end(fwd_hits[1]) + 1 + trim_left
      amp_end <- start(rev_hits[1]) - 1 - trim_right
      if (amp_start < amp_end) {
        amplicon <- subseq(seq, amp_start, amp_end)
        results <- c(results, as.character(amplicon))
      }
    }
  }
  return(unique(results))
}

allo_v3v4 <- extract_v3v4(allobacillus_fasta, fwd_primer, rev_primer, trim_left, trim_right)
imte_v3v4 <- extract_v3v4(imtechella_fasta, fwd_primer, rev_primer, trim_left, trim_right)

if (length(allo_v3v4) == 0) warning("No V3V4 amplicon extracted for Allobacillus")
if (length(imte_v3v4) == 0) warning("No V3V4 amplicon extracted for Imtechella")

message(sprintf("V3V4 refs: Allobacillus %d seqs (%d bp), Imtechella %d seqs (%d bp)",
                length(allo_v3v4), nchar(allo_v3v4[1]),
                length(imte_v3v4), nchar(imte_v3v4[1])))

# --- Fuzzy-match ASVs to spike-in V3V4 references using edit distance ---
asv_seqs <- taxa_names(ps)

min_adist <- function(asvs, refs) {
  if (length(refs) == 0) return(rep(NA_integer_, length(asvs)))
  d <- adist(asvs, refs)
  apply(d, 1, min)
}

dist_allo <- min_adist(asv_seqs, allo_v3v4)
dist_imte <- min_adist(asv_seqs, imte_v3v4)

is_allo <- !is.na(dist_allo) & dist_allo <= max_edit_dist
is_imte <- !is.na(dist_imte) & dist_imte <= max_edit_dist

message(sprintf("Allobacillus: %d ASVs matched (edit dist <= %d)", sum(is_allo), max_edit_dist))
message(sprintf("Imtechella: %d ASVs matched (edit dist <= %d)", sum(is_imte), max_edit_dist))

spike_label <- rep(NA_character_, length(asv_seqs))
spike_label[is_allo] <- "spike_in_Allobacillus_halotolerans"
spike_label[is_imte] <- "spike_in_Imtechella_halotolerans"

# --- Update taxonomy for spike-in ASVs ---
match_spiked <- ps %>% tax_names2rank(colname = "unique")

tt <- as.data.frame(tax_table(match_spiked))
matched_idx <- which(!is.na(spike_label))
tt$Species[matched_idx] <- spike_label[matched_idx]
tax_table(match_spiked) <- tax_table(as.matrix(tt))

# Extract matched spike-in ASVs and apply manual taxonomy
spiked_asvs <- match_spiked %>%
  tax_table() %>%
  as.data.frame() %>%
  filter(Species %in% spiked_in_species) %>%
  select(matched_asv = Species, unique) %>%
  left_join(manual_taxa, by = "matched_asv") %>%
  select(-matched_asv) %>%
  as_tibble()

if (nrow(spiked_asvs) != 0) {
  # Build updated tax_table: non-spike-in taxa + spike-in taxa with manual taxonomy
  updated_tax_table <- match_spiked %>%
    tax_select(spiked_in_species, ranks_searched = "Species", deselect = TRUE) %>%
    tax_table() %>%
    as.data.frame() %>%
    bind_rows(spiked_asvs) %>%
    remove_rownames() %>%
    column_to_rownames("unique") %>%
    as.matrix() %>%
    tax_table()

  tax_table(ps) <- updated_tax_table

  ps %>%
    ps_join(
      ps %>%
        tax_transform(trans = "compositional", rank = "Species") %>%
        ps_otu2samdat(spiked_in_species) %>%
        samdat_tbl() %>%
        rename_with(~ ifelse(str_detect(., "spike_in"), str_c(., "_rel"), .)) %>%
        mutate(sample_id = .sample_name) %>%
        select(sample_id, starts_with("spike_in_"))
    ) %>%
    ps_join(
      ps %>%
        tax_transform(trans = "identity", rank = "Species") %>%
        ps_otu2samdat(spiked_in_species) %>%
        samdat_tbl() %>%
        rename_with(~ ifelse(str_detect(., "spike_in"), str_c(., "_counts"), .)) %>%
        mutate(sample_id = .sample_name) %>%
        select(sample_id, starts_with("spike_in_"))
    ) %>%
    ps_mutate(total_spike_in_rel = rowSums(across(ends_with("_rel")))) %>%
    ps_mutate(total_spike_in_counts = rowSums(across(ends_with("_counts")))) %>%
    ps_mutate(spike_in_ratio = spike_in_Imtechella_halotolerans_counts / spike_in_Allobacillus_halotolerans_counts) %>%
    ps_mutate(distance_from_expected = (3 / 7) - spike_in_ratio) %>%
    tax_select(spiked_in_species, ranks_searched = "Species", deselect = TRUE) %>%
    saveRDS(snakemake@output[["spikein_adjusted_ps"]])
} else {
  # No spike-in ASVs found — save ps with zero spike-in columns
  ps %>%
    ps_mutate(spike_in_Allobacillus_halotolerans_rel = 0) %>%
    ps_mutate(spike_in_Imtechella_halotolerans_rel = 0) %>%
    ps_mutate(spike_in_Allobacillus_halotolerans_counts = 0) %>%
    ps_mutate(spike_in_Imtechella_halotolerans_counts = 0) %>%
    ps_mutate(total_spike_in_rel = 0) %>%
    ps_mutate(total_spike_in_counts = 0) %>%
    saveRDS(snakemake@output[["spikein_adjusted_ps"]])
}
