# Drop ASVs with a long G or C homopolymer run (NextSeq 2-color dark-cycle junk).
#
# In 2-color chemistry a dark cycle (no signal) is base-called G with HIGH quality,
# so neither -q quality trimming nor DADA2's maxEE removes it. A dark forward read
# yields a poly-G ASV; a dark reverse read yields a poly-C ASV after merge/revcomp.
# Real V3V4 16S ASVs cap at ~6 bp G/C homopolymers, so a >=10 bp run cleanly marks
# artifact with no false positives. See investigations/06_polyg_polyc.
#
# Input : chimera-free seqtab (samples x ASVs matrix, colnames = ASV sequences).
# Output: filtered seqtab + a one-row TSV recording what was dropped.

suppressPackageStartupMessages(library(tidyverse))

seqtab <- readRDS(snakemake@input[[1]])
thr    <- as.integer(snakemake@params[["max_run"]])
seqs   <- colnames(seqtab)

# longest run of `base` in each sequence
max_run <- function(s, base) {
  m <- gregexpr(paste0(base, "+"), s)
  vapply(m, function(g) if (g[1] == -1L) 0L else max(attr(g, "match.length")), integer(1))
}

maxG <- max_run(seqs, "G")
maxC <- max_run(seqs, "C")
is_poly <- (maxG >= thr) | (maxC >= thr)

total_reads <- sum(seqtab)
poly_reads  <- sum(seqtab[, is_poly, drop = FALSE])

stats <- tibble(
  threshold        = thr,
  n_asv_total      = length(seqs),
  n_asv_dropped    = sum(is_poly),
  pct_asv_dropped  = round(100 * mean(is_poly), 4),
  reads_total      = total_reads,
  reads_dropped    = poly_reads,
  pct_reads_dropped = round(100 * poly_reads / max(total_reads, 1), 5),
  max_G_run        = max(maxG),
  max_C_run        = max(maxC)
)
write_tsv(stats, snakemake@output[["stats"]])

message(sprintf(
  "homopolymer filter (G/C run >= %d): dropped %d of %d ASVs (%.3f%%), %d of %d reads (%.4f%%); max G-run=%d, max C-run=%d",
  thr, sum(is_poly), length(seqs), 100 * mean(is_poly),
  poly_reads, total_reads, 100 * poly_reads / max(total_reads, 1), max(maxG), max(maxC)))

seqtab_filt <- seqtab[, !is_poly, drop = FALSE]
saveRDS(seqtab_filt, snakemake@output[["seqtab"]])
