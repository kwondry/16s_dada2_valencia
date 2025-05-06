
rule get_speciateit_inputs:
    input:
        spikein_adjusted_ps = "outputs/phyloseq/spikein_adjusted_ps-{run}.rds"
    output:
        asvs = "outputs/cst_processing/speciateit-{run}/asvs.fa",
        count_table = "outputs/cst_processing/speciateit-{run}/count_table.csv"
    conda:
        "../envs/16s_tools.yaml"
    resources:
        cpus_per_task=2, 
        mem_mb=4000,
        runtime="1h",
        partition="short"
    script:
        "../scripts/get_speciateit_input.R"


rule speciateit_classify:
    input:
        asvs = "outputs/cst_processing/speciateit-{run}/asvs.fa",
        speciateit_files = config["speciateit_db"]
    output:
        classified_taxa = "outputs/cst_processing/speciateit-{run}/MC_order7_results.txt"
    params:
        region=lambda wildcards: "vSpeciateIT_V4V4" if run_regions[wildcards.run] == "V4" else "vSpeciateIT_V3V4"
    conda:
        "../envs/16s_tools.yaml"
    resources:
        cpus_per_task=2, 
        mem_mb=4000,
        runtime="1h",
        partition="short"
    shell:
        """
        export PATH="{input.speciateit_files}/bin/linux:$PATH"
        output_dir=$(dirname {output.classified_taxa})
        classify -d {input.speciateit_files}/vSpeciateDB_models/{params.region}   -i {input.asvs} -o $output_dir
        """

rule speciateit_count_table:
    input:
        count_table = "outputs/cst_processing/speciateit-{run}/count_table.csv",
        classified_taxa = "outputs/cst_processing/speciateit-{run}/MC_order7_results.txt",
        speciateit_files = config["speciateit_db"]
    output:
        speciateIT_count_table = "outputs/cst_processing/speciateit-{run}/count_table_speciateIT.csv"
    conda:
        "../envs/16s_tools.yaml"
    resources:
        cpus_per_task=2,
        mem_mb=2000,
        runtime="1h",
        partition="short"
    shell:
        """
        output_dir=$(dirname {output.speciateIT_count_table})
        cd $output_dir
        python {input.speciateit_files}/bin/count_table.py   -c ../../../{input.count_table}   -s ../../../{input.classified_taxa}
        """


rule valencia:
    input:
        speciateIT_count_table = "outputs/cst_processing/speciateit-{run}/count_table_speciateIT.csv",
        speciateit_files = config["speciateit_db"]
    output:
        cst_assignments = "outputs/cst_processing/valencia-{run}/cst_assignments.csv"
    conda:
        "../envs/16s_tools.yaml"
    resources:
        cpus_per_task=2, 
        mem_mb=2000,
        runtime="1h",
        partition="short"
    shell:
        """
        python {input.speciateit_files}/Valencia_v1.1.py   -ref {input.speciateit_files}/VALENCIA2_CST_centroids_19Aug2024.csv   -i {input.speciateIT_count_table}  -o {output.cst_assignments}
        mv {output.cst_assignments}.csv {output.cst_assignments}
        """
        
rule adding_csts_to_ps:
    input:  
        cst_assignments = "outputs/cst_processing/valencia-{run}/cst_assignments.csv",
        spikein_adjusted_ps = "outputs/phyloseq/spikein_adjusted_ps-{run}.rds"
    output:
        ps_with_csts = "outputs/phyloseq/ps_with_cst-{run}.rds"
    conda:
        "../envs/16s_tools.yaml"
    resources:
        cpus_per_task=2, 
        mem_mb=2000,
        runtime="1h",
        partition="short"
    script:
        "../scripts/adding_csts_to_ps.R"


rule get_speciateit_taxa:
    input:
        ps_with_csts = "outputs/phyloseq/ps_with_cst-{run}.rds",
        speciateit_output = "outputs/cst_processing/speciateit-{run}/MC_order7_results.txt",
        speciateit_files = config["speciateit_db"]
    output:
        speciateit_ps = "outputs/results/speciateit_ps-{run}.rds"
    conda:
        "../envs/16s_tools.yaml"
    resources:
        cpus_per_task=2, 
        mem_mb=2000,
        runtime="1h",
        partition="short"
    script:
        "../scripts/get_speciateit_taxa.R"
