library(dada2)
library(tidyverse)

# get all the unique ASVs from this table
asvs <- readRDS(snakemake@input[["seqtab"]]) %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "ASV", values_to = "count") %>%
    group_by(sample) %>%
    ungroup() %>%
    select(ASV) %>%
    distinct() %>%
    pull(ASV)

# run dada2 assign taxonomy and add species
tt <- assignTaxonomy(asvs, "/n/groups/kwon/joseph/dbs/GTDB_bac120_arc53_ssu_r207_fullTaxo.fa.gz", multithread = snakemake@threads, tryRC=TRUE)
tt_plus <- addSpecies(tt, "/n/groups/kwon/joseph/dbs/GTDB_bac120_arc53_ssu_r207_Species.fa.gz", verbose = TRUE, tryRC=TRUE)

# save the gtdb taxonomic calls
saveRDS(tt_plus, snakemake@output[["taxonomy"]])


