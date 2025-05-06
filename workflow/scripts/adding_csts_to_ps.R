#this script adds cst assingments from valencia to a phyloseq's sample metadata

library(tidyverse)
library(phyloseq)
library(microViz)


ps <- readRDS(paste(snakemake@input[["spikein_adjusted_ps"]]))
cst_assignments <- read_csv(paste(snakemake@input[["cst_assignments"]]))

# joins csts
ps_with_csts <- ps_join(
      ps, 
      match_sample_names = "sampleID",
      cst_assignments %>% 
        select(sampleID, CST, subCST, score)
    )

saveRDS(ps_with_csts, snakemake@output[["ps_with_csts"]])
