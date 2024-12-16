# This script adds the spike in relative abundance and number of reads to the sample_data sheet and removes spike in ASVs from the phyloseq.

library(phyloseq)
library(microViz)
library(tidysq)
library(tidyverse)


ps <- readRDS(snakemake@input[["ps"]])

#hardcoded spike in sequences and manually set spike in taxonomy
allobacillus_fasta <- read_fasta("scripts/spike_in_files/D6320.refseq/16S/Allobacillus.halotolerans.16S.fasta") 
imtechella_fasta <- read_fasta("scripts/spike_in_files/D6320.refseq/16S/Imtechella.halotolerans.16S.fasta") 
manual_taxa <- read_csv("scripts/spike_in_files/manual_taxa.csv")
spiked_in_species <- c("spike_in_Allobacillus_halotolerans", "spike_in_Imtechella_halotolerans")


#comparing to spiked_asv
#This renames the species name of an asv to spike_in_.
match_spiked <- ps %>%
                tax_names2rank(colname = "unique") %>%
                tax_mutate(Species = case_when(
                    map_lgl(unique, \(x) any(allobacillus_fasta$sq %has% x)) ~ "spike_in_Allobacillus_halotolerans",
                    map_lgl(unique, \(x) any(imtechella_fasta$sq %has% x)) ~ "spike_in_Imtechella_halotolerans",
                .default = Species))

#asvs and ranks of spiked asv, this renames all ranks of the spiked asv
spiked_asvs <- match_spiked %>%
                tax_table() %>%
                as.data.frame() %>%
                filter(Species %in% spiked_in_species) %>%
                select(matched_asv = Species, unique) %>%
                left_join(manual_taxa, by = "matched_asv") %>%
                select(-matched_asv) %>%
                as_tibble()

#If there are spiked asvs present
if(nrow(spiked_asvs != 0)) {
    #extract a tax table, filter out rows that are the spike in, and then add the 2 from above with udpated all ranks
    updated_tax_table <- match_spiked %>%
                            tax_select(spiked_in_species, ranks_searched = "Species", deselect = TRUE) %>%
                            tax_table() %>%
                            as.data.frame() %>%
                            bind_rows(spiked_asvs) %>%
                            remove_rownames() %>% 
                            column_to_rownames("unique") %>%
                            as.matrix() %>%
                            tax_table()


    #updating ps
    tax_table(ps) <- updated_tax_table

    ps %>%
        ps_join( # adding to samdat the relative abundance of the spike ins
            ps %>%  
                tax_transform(trans = "compositional", rank = "Species") %>%
                ps_otu2samdat(spiked_in_species) %>% 
                samdat_tbl() %>%
                rename_with(~ifelse(str_detect(., "spike_in"), str_c(., "_rel"), .)) %>%
                mutate(sample_id = .sample_name) %>%
                select(sample_id, starts_with("spike_in_"))
        ) %>%
        ps_join( # adding to samdat the spike in counts
            ps %>%
                tax_transform(trans = "identity", rank = "Species") %>%
                ps_otu2samdat(spiked_in_species) %>% 
                samdat_tbl() %>%
                rename_with(~ifelse(str_detect(., "spike_in"), str_c(., "_counts"), .)) %>%
                mutate(sample_id = .sample_name) %>%
                select(sample_id, starts_with("spike_in_"))
        ) %>%
        ps_mutate(total_spike_in_rel = rowSums(across(ends_with("_rel")))) %>%
        ps_mutate(total_spike_in_counts = rowSums(across(ends_with("_counts")))) %>%
        tax_select(spiked_in_species, ranks_searched = "Species", deselect = TRUE) %>%
        saveRDS(snakemake@output[["spikein_adjusted_ps"]])
} else {
    # if no spiked asvs are present, we can leave the ps as is
    ps %>% 
        ps_mutate(spike_in_Allobacillus_halotolerans_rel = 0) %>%
        ps_mutate(spike_in_Imtechella_halotolerans_rel = 0) %>%
        ps_mutate(spike_in_Allobacillus_halotolerans_counts = 0) %>%
        ps_mutate(spike_in_Imtechella_halotolerans_counts = 0) %>%
        ps_mutate(total_spike_in_rel = 0) %>%
        ps_mutate(total_spike_in_counts = 0) %>%
        saveRDS(snakemake@output[["spikein_adjusted_ps"]])

}
