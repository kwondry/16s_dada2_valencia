import pandas as pd
from pathlib import Path
import subprocess

#INPUTS
run_sheet = "sequencing_runsheet.csv"

#path to generated sample sheet
samplesheet_output = "outputs/input/samplesheet.csv"

#creates samplesheet of all samples
if not os.path.exists(samplesheet_output):
    subprocess.run(["python", "scripts/process_samples.py", run_sheet, samplesheet_output], check=True)


#getting pipeline inputs
read_paths = pd.read_csv(samplesheet_output).to_dict(orient="records")
se_read_paths = [rec for rec in read_paths if rec['sequence_type'] == 'single']
pe_read_paths = [rec for rec in read_paths if rec['sequence_type'] == 'paired']
all_runs = set([rec["run_id"] for rec in read_paths])

#getting single and paired runs
se_all_runs = set([rec["run_id"] for rec in se_read_paths])
pe_all_runs = set([rec["run_id"] for rec in pe_read_paths])

sample_map_se = {rec["sample_id"]: rec['fwd'] for rec in read_paths if rec["sequence_type"] == "single"}
sample_map_pe = {rec["sample_id"]: {'fwd': rec['fwd'], 'rev': rec['rev']} for rec in read_paths if rec["sequence_type"] == "paired"}

#gets all samples for single end runs
all_samples_se = {}
for run in se_all_runs: all_samples_se[run] = [rec["sample_id"] for rec in read_paths if rec["run_id"] == run]
se_sample_ids = [s for r in all_samples_se.keys() for s in all_samples_se[r]]

#gets all samples for paired end runs
all_samples_pe = {}
for run in pe_all_runs: all_samples_pe[run] = [rec["sample_id"] for rec in read_paths if rec["run_id"] == run]
pe_sample_ids = [s for r in all_samples_pe.keys() for s in all_samples_pe[r]]

#getting mapping files associated with each run
all_mapping_files = {rec["run_id"]: rec['mapping_file'] for rec in read_paths}

rule target:
    input:
        expand("outputs/dada2_processing/reports/dada2-pe/quality-profile/{run_pe}-{sample_pe}-quality-profile.png", run_pe = pe_all_runs, sample_pe = pe_sample_ids),
        expand("outputs/{run}-multiqc_report.html", run = all_runs)
        # expand("outputs/results/speciateit_ps-{run}.rds", run = all_runs)

include: "rules/processing_se_runs.smk"
include: "rules/processing_pe_runs.smk"


def get_dada2_output(wildcards):
    return "outputs/dada2_processing/results/dada2-se/{}-seqTab.collapsed.RDS".format(wildcards.run) if wildcards.run in se_all_runs else "outputs/dada2_processing/results/dada2-pe/{}-seqTab.collapsed.RDS".format(wildcards.run)

rule gtdb_assign_taxonomy:
    input:
        seqtab = get_dada2_output,
        assign_taxonomy = "/n/groups/kwon/joseph/dbs/GTDB_bac120_arc53_ssu_r207_fullTaxo.fa.gz",
        assign_species = "/n/groups/kwon/joseph/dbs/GTDB_bac120_arc53_ssu_r207_Species.fa.gz"
    output:
        taxonomy = "outputs/dada2_processing/results/dada2/gtdb_taxonomy-{run}.RDS"
    conda:
        "envs/16s_tools.yaml"
    threads:
        8
    resources:
        cpus_per_task=8, 
        mem_mb=16000,
        runtime="8h",
        partition="short"
    script:
        "scripts/assign_taxonomy.R"

        
def get_mapping(wildcards):
    return all_mapping_files[wildcards.run]

rule making_phyloseq:
    input:
        mapping_file= get_mapping,
        seqtab = get_dada2_output,
        taxonomy = "outputs/dada2_processing/results/dada2/gtdb_taxonomy-{run}.RDS"
    output:
        ps = "outputs/results/gtdb_ps-{run}.rds"
    params:
        threshold = 0
    conda:
        "envs/16s_tools.yaml"
    resources:
        cpus_per_task=2, 
        mem_mb=4000,
        runtime="8h",
        partition="short"
    script:
        "scripts/making_phyloseq.R"

# not on conda
rule install_tidysq:
    output:
        tidy_install = touch("outputs/flags/tidyq_installed.done")
    conda:
        "envs/16s_tools.yaml"
    resources:
        cpus_per_task=2, 
        mem_mb=1000,
        runtime="2h",
        partition="short"
    shell:
        """
        Rscript -e 'if(!require(tidysq)){{ install.packages("tidysq", repos="https://cloud.r-project.org", quiet = TRUE) }}'
        echo "tidysq successfully installed" > {output.tidy_install}
        """

rule spikein_adjustment:
    input:
        tidy_install = "outputs/flags/tidyq_installed.done",
        ps = "outputs/results/gtdb_ps-{run}.rds"
    output: 
        spikein_adjusted_ps = "outputs/phyloseq/spikein_adjusted_ps-{run}.rds"
    conda:
        "envs/16s_tools.yaml"
    resources:
        cpus_per_task=2, 
        mem_mb=8000,
        runtime="8h",
        partition="short"
    script:
        "scripts/spikein_ps_adjustment.R"


include: "rules/speciateit_valencia.smk"


rule get_multiqc_inputs:
    input:
        speciateit_ps = "outputs/results/speciateit_ps-{run}.rds"
    output:
        filt_read_counts = "outputs/results/{run}-read_counts.csv",
        ordplot = "outputs/results/{run}-asv_ordination.png",
        barplot = "outputs/results/{run}-abundance_barplots.png",
    conda:
        "envs/16s_tools.yaml"
    params:
        run = lambda wildcards: wildcards.run
    resources:
        cpus_per_task=2, 
        mem_mb=8000,
        runtime="8h",
        partition="short"
    script:
        "scripts/get_multiqc_inputs.R"

# generates one per run
rule multiqc:
    input:
        runs = "outputs/results/speciateit_ps-{run}.rds",
        ordplot = "outputs/results/{run}-asv_ordination.png",
        barplot = "outputs/results/{run}-abundance_barplots.png",
        filt_read_counts = "outputs/results/{run}-read_counts.csv",
        hisat = "outputs"
    output:
        report="outputs/{run}-multiqc_report.html"
    params:
        use_input_files_only=True,
        extra="-c configs/mqc_config.yaml"
    resources:
        cpus_per_task=2, 
        mem_mb=8000,
        runtime="8h",
        partition="short"
    log:
        "outputs/logs/{run}-multiqc.log",
    wrapper:
        "v5.2.1/bio/multiqc"
