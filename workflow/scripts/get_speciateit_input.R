# this script generates inputs for speciateit
library(microViz)
library(tidyverse)

ps <- readRDS(snakemake@input[["spikein_adjusted_ps"]])


asvs_to_write <- tibble(asv=microViz::otu_get(ps) %>% 
                  colnames()) %>% 
                  mutate(fasta_id=str_c(">",asv, "\n")) %>% 
                  mutate(asv = str_c(asv, "\n"))

asvs_to_write %>%
  mutate(line_to_write=map2_chr(fasta_id, asv, str_c)) %>%
  pull(line_to_write) %>%
  reduce(str_c) %>%
  write_file(snakemake@output[["asvs"]])

ps %>%
  otu_get() %>%
  as.data.frame() %>%
  write.csv(snakemake@output[["count_table"]]) # using base R write.csv since we need the rownames