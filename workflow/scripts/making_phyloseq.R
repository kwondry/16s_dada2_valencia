
# This is a script that takes dada2 output, appends the sequencing run's mapping file as the sample data, and transposes the otu_table, so the ASVs are column headers for following script and speciateIt compatibility.
library(phyloseq)
library(microViz)
library(tidyverse)

# Read the ID mapping file
id_mapping <- read_tsv(snakemake@input[["id_mapping"]])

dada2_taxonomy <- readRDS(snakemake@input[["taxonomy"]]) %>%
    as.data.frame() %>%
    as_tibble(rownames = "ASV", .name_repair = "unique") %>%
    select(-c("Species...8")) %>%
    rename(Domain = "Kingdom", Species ="Species...7") %>%
    mutate(across(-ASV, ~ if_else(is.na(.), "unassigned", .)))


count_matrix <- readRDS(snakemake@input[["seqtab"]]) %>%
    as.data.frame() %>%
    rownames_to_column("short_id") %>%
    # Map short IDs back to original IDs
    left_join(id_mapping) %>%
    select(-short_id) %>%
    pivot_longer(-original_id, names_to = "ASV", values_to = "count") %>%
    group_by(original_id) %>%
    filter(sum(count) > snakemake@params[["threshold"]]) %>% # filter out samples with < reads 
    ungroup() %>%
    filter(count > 0) %>%
    pivot_wider(names_from = ASV, values_from = count, values_fill = 0) %>%
    mutate(original_id = str_replace_all(original_id, "_S\\d+(_L00\\d)?_R1_00\\d", "")) %>%
    column_to_rownames("original_id") %>%
    as.matrix()

samp_data <- read_delim(snakemake@input[["mapping_file"]]) %>%
    mutate(`#SampleID`=str_replace_all(`#SampleID`, "\\.", "_")) %>%
    filter(`#SampleID` %in% rownames(count_matrix)) %>%  # filter out samples with set read threshold 
    column_to_rownames("#SampleID") %>%
    mutate(sample_id = rownames(.))
    

# this filters the taxonomy to only ASVs that are in the seqtab
taxonomy_matrix <- count_matrix %>%
    as_tibble(rownames = "sample") %>% 
    pivot_longer(-sample, names_to = "ASV", values_to = "count") %>%
    select(ASV) %>%
    distinct() %>%
    inner_join(dada2_taxonomy) %>%
    select(ASV, Domain, Phylum, Class, Order, Family, Genus, Species) %>%
    column_to_rownames("ASV") %>%
    as.matrix()

    
ps <- phyloseq(otu_table(count_matrix, taxa_are_rows = FALSE), tax_table(taxonomy_matrix), sample_data(samp_data)) %>%
        tax_fix(unknowns = c("unassigned"), anon_unique = FALSE) %>%
        tax_mutate(Genus_Species = str_c(Genus, Species, sep = " "))


saveRDS(ps, snakemake@output[["ps"]])
