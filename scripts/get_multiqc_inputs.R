# this script generates QC plots and files for the multiqc report

library(tidyverse)
library(phyloseq)
library(microViz)
library(svglite)
library(cowplot)

#inputs
ps <- readRDS(paste(snakemake@input[["speciateit_ps"]]))

#generates read count list
list.files("outputs/dada2_processing/reports/dada2-pe/filter-trim-pe", recursive = TRUE, pattern = "\\.tsv$", full.names = TRUE) %>%
  map_df(read_tsv) %>%
  as_tibble() %>%
  filter(str_detect(reads.in, snakemake@params[["run"]])) %>%
  rename(sample_id = reads.in, reads = reads.out) %>%
  separate_wider_delim(reads, names = c("reads_in", "Reads post-filtering"), delim = "\t") %>%
  select(sample_id, `Reads post-filtering`) %>%
  arrange(`Reads post-filtering`) %>%
  mutate(sample_id = str_replace(sample_id, paste0(snakemake@params[["run"]], "-"), "")) %>%
  mutate(sample_id = str_replace(sample_id, ".1.fastq.gz", "")) %>%
  mutate(sample_id = paste0('"',sample_id, '"')) %>%
  write.csv(paste0(snakemake@output["filt_read_counts"]), quote =FALSE, row.names = FALSE)


# ord plot
ordplot <- ps %>%
  ps_mutate(Researcher_Description = str_c(Researcher, Description, sep = "_")) %>%
  tax_transform(rank = "unique", trans = "compositional") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(method = "auto") %>% 
  ord_plot(
  axes = c(1, 2),
  color = "subCST", shape = "Researcher_Description",
  size = 0.05
  ) +
  theme_cowplot(3)+
  background_grid()

#bar plot
barplot <- ps %>%
  ps_mutate(Researcher_Description = str_c(Researcher, Description, sep = "_")) %>%
  comp_barplot(
    tax_level = "Species",
    label = NULL, # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    facet_by = "Researcher_Description",
    bar_width = 0.9,
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) +
  coord_flip() +
  theme_cowplot(3) +
  background_grid()


save_plot(snakemake@output[["ordplot"]], ordplot, base_width = 3, base_height = 2)
save_plot(snakemake@output[["barplot"]], barplot, base_width = 3, base_height = 2.5)

