
def get_reads_pe(wildcards):
    return {"fwd": sample_map_pe[wildcards.sample_pe]['fwd'], "rev": sample_map_pe[wildcards.sample_pe].get('rev', None)}

#assuming all pe runs are done on v3v4 region
rule trim_primers:
    input:
        unpack(get_reads_pe)
    output:
        forward_trimmed="outputs/v3v4_trimmed/{run_pe}-{sample_pe}_R1.trimmed.fastq.gz",
        reverse_trimmed="outputs/v3v4_trimmed/{run_pe}-{sample_pe}_R2.trimmed.fastq.gz"
    params:
        forward_primer=config["primers"]["v3v4"]["forward"],
        reverse_primer=config["primers"]["v3v4"]["reverse"]
    log:
        "outputs/logs/cutadapt_{run_pe}/{run_pe}-{sample_pe}_cutadapt.log"
    conda:
        "../envs/cutadapt.yaml"
    resources:
        cpus_per_task=1, 
        mem_mb=1000,
        runtime="1h",
        partition="short"
    shell:
        """
        cutadapt \
            -g {params.forward_primer} \
            -G {params.reverse_primer} \
            -o {output.forward_trimmed} \
            -p {output.reverse_trimmed} \
            --discard-untrimmed \
            {input.fwd} {input.rev} \
            > {log} 2>&1
        """

rule hisat2_remove_human_sensitive_pe:
    input:
        forward_trimmed="outputs/v3v4_trimmed/{run_pe}-{sample_pe}_R1.trimmed.fastq.gz",
        reverse_trimmed="outputs/v3v4_trimmed/{run_pe}-{sample_pe}_R2.trimmed.fastq.gz"
        # idx_files="/n/groups/kwon/joseph/dbs/combined_T2T_CRCh38_reference_for_host_filtering_hisat2"
    output:
        fwd= "outputs/human_filtered-pe/{run_pe}-{sample_pe}.1.fastq.gz",
        rev= "outputs/human_filtered-pe/{run_pe}-{sample_pe}.2.fastq.gz",
        report="outputs/human_filtered-pe/{run_pe}_reports/{sample_pe}.txt"
    params:
        hisat2_index=config["hisat_db"]
    threads: 
        4
    conda:
        "../envs/hisat2.yaml"
    resources:
        cpus_per_task=1, 
        mem_mb=16000,
        runtime="1h",
        partition="short"
    shell:
        """
        hisat2 -x {params.hisat2_index} -1 {input.forward_trimmed} -2 {input.reverse_trimmed} \
            -p {threads} --no-spliced-alignment \
            --score-min L,0,-0.6 \
            --un-conc-gz outputs/human_filtered-pe/{wildcards.run_pe}-{wildcards.sample_pe}.%.fastq.gz \
            -S /dev/null 2> {output.report}
        """

rule dada2_quality_profile_pe:
    input:
        # FASTQ file without primer sequences
        expand("outputs/human_filtered-pe/{{run_pe}}-{{sample_pe}}.{orientation}.fastq.gz",orientation=[1,2])
    output:
        "outputs/dada2_processing/reports/dada2-pe/quality-profile/{run_pe}-{sample_pe}-quality-profile.png"
    log:
        "outputs/dada2_processing/logs/dada2-pe/quality-profile/{run_pe}-{sample_pe}-quality-profile-pe.log"
    resources:
        cpus_per_task=1, 
        mem_mb=4000,
        runtime="1h",
        partition="short"
    wrapper:
        "v5.2.1/bio/dada2/quality-profile/wrapper.R"

rule dada2_filter_trim_pe:
    input:
        fwd= "outputs/human_filtered-pe/{run_pe}-{sample_pe}.1.fastq.gz",
        rev= "outputs/human_filtered-pe/{run_pe}-{sample_pe}.2.fastq.gz"
    output:
        filt = "outputs/dada2_processing/filtered-pe/{run_pe}-{sample_pe}.1.fastq.gz",
        filt_rev = "outputs/dada2_processing/filtered-pe/{run_pe}-{sample_pe}.2.fastq.gz",
        stats = "outputs/dada2_processing/reports/dada2-pe/filter-trim-pe/{run_pe}-{sample_pe}.tsv"
    params:
        maxEE=config["dada2"]["pe"]["maxEE"],
        truncLen=config["dada2"]["pe"]["truncLen"],
        trimLeft=config["dada2"]["pe"]["trimLeft"]
    log:
        "outputs/logs/dada2-pe/filter-trim-pe/{run_pe}-{sample_pe}.log"
    threads: 
        4
    resources:
        cpus_per_task=4, 
        mem_mb=4000,
        runtime="1h",
        partition="short"
    wrapper:
       "v5.2.1/bio/dada2/filter-trim/wrapper.R"

        
rule dada2_learn_errors_pe:
    input:
        lambda wildcards: expand("outputs/dada2_processing/filtered-pe/{{run_pe}}-{sample_pe}.{{orientation}}.fastq.gz", sample_pe=all_samples_pe[wildcards.run_pe])
    output:
        err = "outputs/dada2_processing/results/dada2-pe/model_{run_pe}_{orientation}.RDS",# save the error model
        plot = "outputs/dada2_processing/reports/dada2-pe/errors_{run_pe}_{orientation}.png",# plot observed and estimated rates
    params:
        randomize = True
    log:
        "outputs/logs/dada2-pe/learn-errors/learn-errors_{run_pe}_{orientation}.log"
    threads: 
        4
    resources:
        cpus_per_task=4,
        mem_mb=4000,
        runtime="8h",
        partition="short"
    wrapper:
        "v5.2.1/bio/dada2/learn-errors/wrapper.R"


rule dada2_dereplicate_fastq_pe:
    input:
        "outputs/dada2_processing/filtered-pe/{run_pe}-{sample_pe}.{orientation}.fastq.gz"
    output:
        "outputs/dada2_processing/uniques-pe/{run_pe}-{sample_pe}.{orientation}.RDS"
    log:
        "outputs/logs/dada2/dereplicate-fastq/{run_pe}-{sample_pe}.{orientation}.log"
    resources:
        cpus_per_task=1, 
        mem_mb=4000,
        runtime="1h",
        partition="short"
    wrapper:
        "v5.2.1/bio/dada2/dereplicate-fastq/wrapper.R"


rule dada2_sample_inference_pe:
    input:
        derep="outputs/dada2_processing/uniques-pe/{run_pe}-{sample_pe}.{orientation}.RDS",
        err="outputs/dada2_processing/results/dada2-pe/model_{run_pe}_{orientation}.RDS" # Error model
    output:
        "outputs/dada2_processing/denoised-pe/{run_pe}-{sample_pe}.{orientation}.RDS" # Inferred sample composition
    log:
        "outputs/logs/dada2/sample-inference/{run_pe}-{sample_pe}.{orientation}.log"
    threads: 
        4
    resources:
        cpus_per_task=4, 
        mem_mb=8000,
        runtime="1h",
        partition="short"
    wrapper:
        "v5.2.1/bio/dada2/sample-inference/wrapper.R"

rule dada2_merge_pairs:
    input:
      dadaF="outputs/dada2_processing/denoised-pe/{run_pe}-{sample_pe}.1.RDS",# Inferred composition
      dadaR="outputs/dada2_processing/denoised-pe/{run_pe}-{sample_pe}.2.RDS",
      derepF="outputs/dada2_processing/uniques-pe/{run_pe}-{sample_pe}.1.RDS",# Dereplicated sequences
      derepR="outputs/dada2_processing/uniques-pe/{run_pe}-{sample_pe}.2.RDS"
    output:
        "outputs/dada2_processing/merged/{run_pe}-{sample_pe}.RDS"
    log:
        "outputs/logs/dada2/merge-pairs/{run_pe}-{sample_pe}.log"
    threads:
        4
    params:
        minOverlap=config["dada2"]["pe"]["merge"]["minOverlap"],
        maxMismatch=config["dada2"]["pe"]["merge"]["maxMismatch"],
    resources:
        cpus_per_task=4, 
        mem_mb=8000,
        runtime="1h",
        partition="short"
    wrapper:
        "v5.2.1/bio/dada2/merge-pairs/wrapper.R"

rule dada2_make_table_pe:
    input:
        lambda wildcards: expand("outputs/dada2_processing/merged/{{run_pe}}-{sample_pe}.RDS", sample_pe=all_samples_pe[wildcards.run_pe])
    output:
        "outputs/dada2_processing/results/dada2-pe/{run_pe}-seqTab-pe.RDS"
    params:
        names = lambda wildcards: all_samples_pe[wildcards.run_pe]
    log:
        "outputs/logs/dada2/make-table/{run_pe}-make-table-pe.log"
    threads:
        1
    resources:
        cpus_per_task=1, 
        mem_mb=8000,
        runtime="1h",
        partition="short"
    wrapper:
        "v5.2.1/bio/dada2/make-table/wrapper.R"

rule dada2_remove_chimeras_pe:
    input:
        "outputs/dada2_processing/results/dada2-pe/{run_pe}-seqTab-pe.RDS" # Sequence table
    output:
        "outputs/dada2_processing/results/dada2-pe/{run_pe}-seqTab.nochimeras.RDS" # Chimera-free sequence table
    log:
        "outputs/logs/dada2/remove-chimeras/{run_pe}-remove-chimeras.log"
    threads:
        16
    resources:
        cpus_per_task=16, 
        mem_mb=8000,
        runtime="24h",
        partition="medium"
    wrapper:
        "v5.2.1/bio/dada2/remove-chimeras/wrapper.R"

# rule dada2_collapse_nomismatch_pe:
#     input:
#         "outputs/dada2_processing/results/dada2-pe/{run_pe}-seqTab.nochimeras.RDS" # Chimera-free sequence table
#     output:
#         "outputs/dada2_processing/results/dada2-pe/{run_pe}-seqTab.collapsed.RDS"
#     log:
#         "outputs/logs/dada2/collapse-nomismatch/{run_pe}-collapse-nomismatch.log"
#     threads:
#         1
#     resources:
#         cpus_per_task=1, 
#         mem_mb=8000,
#         runtime="12h",
#         partition="short"
#     wrapper:
#         "v5.2.1/bio/dada2/collapse-nomismatch/wrapper.R"
