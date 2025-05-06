# this script generates QC plots and files for the multiqc report

library(tidyverse)
library(phyloseq)
library(microViz)
library(svglite)
library(cowplot)

#inputs
ps <- readRDS(paste(snakemake@input[["speciateit_ps"]]))

#generates read count list
list.files(snakemake@params[["read_counts"]], recursive = TRUE, pattern = "\\.tsv$", full.names = TRUE) %>% 
  as_tibble() %>%
  filter(str_detect(value, snakemake@params[["run"]])) %>%
  map_df(read_tsv, col_types = cols(.default = col_character())) %>%
  rename(sample_id = reads.in, reads = reads.out) %>%
  separate_wider_delim(reads, names = c("Reads post cutadapt and human filtering", "Reads post dada2 filter-trim"), delim = "\t") %>%
  select(sample_id, `Reads post cutadapt and human filtering`,`Reads post dada2 filter-trim`) %>%
  bind_rows(
    read_delim(snakemake@params[["low_read_samples"]], col_names = c("run","sample_id")) %>%
      filter(str_detect(run, snakemake@params[["run"]])) %>%
      select(-run) %>%
      mutate(sample_id = str_replace(sample_id, "_R(1|2)_001.fastq", ""),
            `Reads post cutadapt and human filtering` = "<250",
            `Reads post dada2 filter-trim` = "-") %>%
      mutate(sample_id = str_replace(sample_id, ".*/", ""))
  ) %>%
  arrange(`Reads post dada2 filter-trim`) %>%
  mutate(sample_id = str_replace(sample_id, paste0(snakemake@params[["run"]], "-"), "")) %>%
  mutate(sample_id = str_replace(sample_id, ".1.fastq.gz", "")) %>%
  mutate(sample_id = paste0('"',sample_id, '"')) %>%
  write.csv(paste0(snakemake@output["filt_read_counts"]), quote =FALSE, row.names = FALSE)


# ord plot
ordplot <- ps %>%
  #ps_mutate(Researcher_Description = str_c(Researcher, Description, sep = "_")) %>%
  tax_transform(rank = "unique", trans = "compositional") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(method = "auto") %>% 
  ord_plot(
  color = "subCST",
  size = 0.05
  ) +
  theme_cowplot(3)+
  background_grid()

#bar plot
barplot <- ps %>%
  # ps_mutate(Researcher_Description = str_c(Researcher, Description, sep = "_")) %>%
  comp_barplot(
    tax_level = "Species",
    label = NULL, # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    #facet_by = "Researcher_Description",
    bar_width = 0.9,
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) +
  coord_flip() +
  theme_cowplot(3) +
  background_grid()


save_plot(snakemake@output[["ordplot"]], ordplot, base_width = 2.5, base_height = 1.5, dpi = 300)
save_plot(snakemake@output[["barplot"]], barplot, base_width = 3, base_height = 2.5, dpi = 300)

