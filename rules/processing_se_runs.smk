
def get_reads_se(wildcards):
    return {"forward_reads":sample_map_se[wildcards.sample_se]}

rule hisat2_remove_human_sensitive_se:
    input:
        unpack(get_reads_se)
    output:
        fwd= "outputs/human_filtered-se/{run_se}-{sample_se}.1.fastq.gz",
        report="outputs/human_filtered-se/reports/{run_se}-{sample_se}_hisat2_report.txt"
    params:
        hisat2_index=config["hisat_db"]
    threads:
        4
    conda:
        "../envs/hisat2.yaml"
    resources:
        cpus_per_task=4, 
        mem_mb=16000,
        runtime="1h",
        partition="short"
    shell:
        """
        hisat2 -x {params.hisat2_index} -U {input.forward_reads} \
            -p {threads} --no-spliced-alignment \
            --score-min L,0,-0.6 \
            --un-gz outputs/human_filtered-se/{wildcards.run_se}-{wildcards.sample_se}.1.fastq.gz \
            -S /dev/null 2> {output.report}
        """

rule dada2_filter_trim_se:
    input:
        fwd= "outputs/human_filtered-se/{run_se}-{sample_se}.1.fastq.gz",
    output:
        filt = "outputs/dada2_processing/filtered-se/{run_se}-{sample_se}.1.fastq.gz",
        stats = "outputs/dada2_processing/reports/dada2-se/filter-trim-se/{run_se}-{sample_se}.tsv"
    params:
        maxEE=1,
        truncLen=230,
        trimLeft=10
    log:
        "outputs/logs/dada2-se/filter-trim-se/{run_se}-{sample_se}.log"
    threads: 
        4
    resources:
        cpus_per_task=4, 
        mem_mb=4000,
        runtime="1h",
        partition="short"
    wrapper:
       "v5.2.1/bio/dada2/filter-trim/wrapper.R"


rule dada2_learn_errors_se:
    input:
        lambda wildcards: expand("outputs/dada2_processing/filtered-se/{{run_se}}-{sample_se}.1.fastq.gz", sample_se=all_samples_se[wildcards.run_se])
    output:
        err = "outputs/dada2_processing/results/dada2-se/model_{run_se}.RDS",# save the error model
        plot = "outputs/dada2_processing/reports/dada2-se/errors_{run_se}.png",# plot observed and estimated rates
    params:
        randomize = True
    log:
        "outputs/logs/dada2-se/learn-errors/learn-errors_{run_se}.log"
    threads: 
        4
    resources:
        cpus_per_task=4, 
        mem_mb=4000,
        runtime="1h",
        partition="short"
    wrapper:
        "v5.2.1/bio/dada2/learn-errors/wrapper.R"


rule dada2_dereplicate_fastq_se:
    input:
        "outputs/dada2_processing/filtered-se/{run_se}-{sample_se}.1.fastq.gz"
    output:
        "outputs/dada2_processing/uniques-se/{run_se}-{sample_se}.RDS"
    log:
        "outputs/logs/dada2-se/dereplicate-fastq/{run_se}-{sample_se}.log"
    threads: 
        4
    resources:
        cpus_per_task=4, 
        mem_mb=4000,
        runtime="1h",
        partition="short"
    wrapper:
        "v5.2.1/bio/dada2/dereplicate-fastq/wrapper.R"

rule dada2_sample_inference_se:
    input:
        derep = "outputs/dada2_processing/uniques-se/{run_se}-{sample_se}.RDS",
        err = "outputs/dada2_processing/results/dada2-se/model_{run_se}.RDS" # Error model
    output:
        "outputs/dada2_processing/denoised-se/{run_se}-{sample_se}.1.RDS" # Inferred sample composition
    log:
        "outputs/logs/dada2-se/sample-inference/{run_se}-{sample_se}.log"
    threads: 
        4
    resources:
        cpus_per_task=4, 
        mem_mb=4000,
        runtime="1h",
        partition="short"
    wrapper:
        "v5.2.1/bio/dada2/sample-inference/wrapper.R"

rule dada2_make_table_se:
    input:
        lambda wildcards: expand("outputs/dada2_processing/denoised-se/{{run_se}}-{sample_se}.1.RDS", sample_se=all_samples_se[wildcards.run_se])
    output:
        "outputs/dada2_processing/results/dada2-se/{run_se}-seqTab-se.RDS"
    log:
        "outputs/logs/dada2-se/make-table/{run_se}-make-table-se.log"
    params:
        names= lambda wildcards: [all_samples_se[wildcards.run_se]] # Sample names instead of paths
    threads: 
        4
    resources:
        cpus_per_task=4, 
        mem_mb=4000,
        runtime="1h",
        partition="short"
    wrapper:
        "v5.2.1/bio/dada2/make-table/wrapper.R"


rule dada2_remove_chimeras_se:
    input:
        "outputs/dada2_processing/results/dada2-se/{run_se}-seqTab-se.RDS" # Sequence table
    output:
        "outputs/dada2_processing/results/dada2-se/{run_se}-seqTab.nochimeras-se.RDS" # Chimera-free sequence table
    log:
        "outputs/logs/dada2-se/remove-chimeras/{run_se}-remove-chimeras-.log"
    threads:
        4
    resources:
        cpus_per_task=4, 
        mem_mb=4000,
        runtime="3h",
        partition="short"
    wrapper:
        "v5.2.1/bio/dada2/remove-chimeras/wrapper.R"

rule dada2_collapse_nomismatch_se:
    input:
        "outputs/dada2_processing/results/dada2-se/{run_se}-seqTab.nochimeras-se.RDS" # Chimera-free sequence table
    output:
        "outputs/dada2_processing/results/dada2-se/{run_se}-seqTab.collapsed.RDS"
    log:
        "outputs/logs/dada2-se/collapse-nomismatch/{run_se}-collapse-nomismatch-.log"
    threads:
        4
    resources:
        cpus_per_task=4, 
        mem_mb=4000,
        runtime="3h",
        partition="short"
    wrapper:
        "v5.2.1/bio/dada2/collapse-nomismatch/wrapper.R"
