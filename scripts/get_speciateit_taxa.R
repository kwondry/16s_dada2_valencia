library(tidyverse)
library(phyloseq)
library(microViz)

#inputs
ps <- readRDS(paste(snakemake@input[["ps_with_csts"]]))
path_to_speciateit_output <- paste(snakemake@input[["speciateit_output"]])
path_to_speciateitfiles <- paste(snakemake@input[["speciateit_files"]])

# function that appends speciateit taxa to a phyloseq
assigning_speciateIT_taxa <- function(ps_to_add_taxa, path_to_speciateit_output) {
  sp_db_lineage <- read_delim(str_c(path_to_speciateitfiles, "/training_data/vSpeciateIT_V4V4.lineage"), col_names = FALSE) %>%
    rbind(
      read_delim(str_c(path_to_speciateitfiles, "/training_data/vSpeciateIT_V3V4.lineage"), col_names = FALSE)
    ) %>%
    rbind(
      read_delim(str_c(path_to_speciateitfiles, "/training_data/vSpeciateIT_V1V3.lineage"), col_names = FALSE)
    ) %>%
    rbind(
    read_delim(str_c(path_to_speciateitfiles, "/training_data/extended_lineages.lineage"), col_names = FALSE)
    ) %>%
    select(Domain = X7, Phylum = X6,  Class = X5, Order = X4, Family = X3, Genus = X2, Species = X1) %>%
    distinct() %>%
    mutate(taxa_id = str_c("taxa_", row_number()))

  #aggregate taxa are specific ranks
  spv4_db_agg <- sp_db_lineage %>%
    mutate(
      agg_domain = Domain,
      agg_phylum = str_c(Domain, Phylum, sep = "~"),
      agg_class = str_c(Domain, Phylum, Class, sep = "~"),
      agg_order = str_c(Domain, Phylum, Class, Order, sep = "~"),
      agg_family = str_c(Domain, Phylum, Class, Order, Family, sep = "~"),
      agg_genus = str_c(Domain, Phylum, Class, Order, Family, Genus,  sep = "~"),
      agg_species = str_c(Domain, Phylum, Class, Order, Family, Genus, Species, sep = "~")
    )


  spv4_db_longer <- sp_db_lineage %>%
    pivot_longer(-taxa_id, names_to = "rank", values_to = "taxa") %>%
    left_join(spv4_db_agg, by = "taxa_id") %>%
    select(rank, taxa, starts_with("agg"))

  ### here is the actual specitateit output and I want to add the full rank
  sp_output_taxonomy <- read_delim(path_to_speciateit_output, col_names = FALSE) %>%
    select(asv = X1, taxa_rank = X2, score = X3, count = X4) %>%
    left_join(
      spv4_db_longer, by = c("taxa_rank" = "taxa"), multiple = "first"
    ) %>% 
    mutate(updated_rank = case_when(
      rank == "Domain" ~ agg_domain,
      rank == "Phylum" ~ agg_phylum,
      rank == "Class" ~ agg_class,
      rank == "Order" ~ agg_order,
      rank == "Family" ~ agg_family,
      rank == "Genus" ~ agg_genus,
      rank == "Species" ~ agg_species
    )) %>%
    select(-starts_with("agg"),-taxa_rank, -score, -count, -rank) %>%
    separate(updated_rank, c("Domain", "Phylum", "Class", "Order", "Family","Genus","Species"), sep = "~") %>%
    column_to_rownames("asv") %>%
    as.matrix()

tax_table(ps_to_add_taxa) <- tax_table(sp_output_taxonomy)

ps_to_add_taxa %>% 
  tax_fix()
}


#changing taxonomy of ps
ps_with_csts_speciateIT_taxa <- ps %>% 
    assigning_speciateIT_taxa(., path_to_speciateit_output)

saveRDS(ps_with_csts_speciateIT_taxa, paste(snakemake@output[["speciateit_ps"]]))