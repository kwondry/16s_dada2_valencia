#!/usr/bin/env Rscript
# Reference-based chimera REVIEW (flag, do NOT delete) for one run.
#
# DADA2's removeBimeraDenovo is conservative and gets weak when the parent pool
# is large (our ASV inflation makes it large). This runs vsearch --uchime_ref
# against the GTDB SSU FASTA to catch chimeras that survive, and writes a review
# table so you can confirm no real vaginal taxa are being flagged BEFORE wiring
# deletion into the DAG. See investigations/05_high_abund_bad_taxonomy and
# investigations/06_polyg_polyc.
#
# This is intentionally standalone (not a default rule): run it on one run,
# inspect <prefix>-chimera_review.tsv (especially the abundant flagged ASVs),
# then we promote it to a deletion rule mirroring filter_homopolymer.R.
#
# Usage:
#   Rscript chimera_ref_review.R <seqTab.filtered.RDS> <gtdb_ssu.fa.gz> <out_prefix> [threads]
# Requires vsearch on PATH (conda: workflow/envs/vsearch.yaml).

suppressPackageStartupMessages({ library(tidyverse) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: chimera_ref_review.R <seqTab.RDS> <gtdb.fa.gz> <out_prefix> [threads]")
seqtab_path <- args[[1]]; db <- args[[2]]; prefix <- args[[3]]
threads <- if (length(args) >= 4) args[[4]] else "4"

seqtab <- readRDS(seqtab_path)
seqs  <- colnames(seqtab)
ids   <- paste0("ASV", seq_along(seqs))
reads <- colSums(seqtab)
prev  <- colSums(seqtab > 0)

# vsearch uses ;size= for the uchime_ref abundance ordering. Write FASTA with
# base R (no Biostrings) so this runs with just tidyverse + vsearch on PATH.
fasta <- paste0(prefix, "-asvs.fasta")
writeLines(paste0(">", ids, ";size=", reads, "\n", seqs), fasta)

uchimeout <- paste0(prefix, "-uchime_ref.txt")
chim_fa   <- paste0(prefix, "-chimeras.fasta")
cmd <- sprintf(
  "vsearch --uchime_ref %s --db %s --uchimeout %s --chimeras %s --threads %s",
  shQuote(fasta), shQuote(db), shQuote(uchimeout), shQuote(chim_fa), threads)
message("Running: ", cmd)
if (system(cmd) != 0) stop("vsearch failed")

# UCHIME tab format: col1 = score, col2 = query (label;size=N), last col = Y/N/? verdict
uc_raw <- read_tsv(uchimeout, col_names = FALSE, show_col_types = FALSE)
last <- ncol(uc_raw)
uc <- tibble(
  id      = str_remove(uc_raw$X2, ";size=.*"),
  score   = uc_raw$X1,
  verdict = uc_raw[[last]])

tbl <- tibble(id = ids, len = nchar(seqs), reads = reads, prevalence = prev, seq = seqs) %>%
  left_join(uc, by = "id") %>%
  mutate(verdict = replace_na(verdict, "N")) %>%
  arrange(desc(verdict == "Y"), desc(reads))

write_tsv(tbl, paste0(prefix, "-chimera_review.tsv"))

flagged <- tbl %>% filter(verdict == "Y")
cat(sprintf("\nFlagged chimeric: %d of %d ASVs (%.2f%%), %s of %s reads (%.3f%%)\n",
            nrow(flagged), nrow(tbl), 100 * nrow(flagged) / nrow(tbl),
            format(sum(flagged$reads), big.mark=","), format(sum(tbl$reads), big.mark=","),
            100 * sum(flagged$reads) / sum(tbl$reads)))
cat("\nMost ABUNDANT flagged ASVs — eyeball these for real-taxon false positives:\n")
print(flagged %>% slice_max(reads, n = 25) %>% select(id, len, reads, prevalence, score), n = 25)
cat(sprintf("\nFull table: %s-chimera_review.tsv\n", prefix))
