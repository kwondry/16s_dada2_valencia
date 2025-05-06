#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=4GB
#SBATCH --partition=priority
#SBATCH --job-name=submit_jobs


# snakemake --executor slurm --default-resources slurm_partition=short --set-resources --configfile configs/config.yaml --use-conda --conda-prefix=~/snakemake_conda_prefix --jobname {rulename}.{jobid} --jobs 2000 --keep-going --rerun-incomplete
snakemake --executor slurm --default-resources slurm_partition=short --set-resources dada2_remove_chimeras_pe:slurm_partition=medium --configfile configs/config.yaml --use-conda --conda-prefix=~/snakemake_conda_prefix --jobname {rulename}.{jobid} --jobs 6000 --keep-going --rerun-incomplete \
--groups trim_primers=group_trim_primers hisat2_remove_human_sensitive_pe=group_hisat2_remove_human_sensitive_pe dada2_quality_profile_pe=group_dada2_quality_profile_pe dada2_filter_trim_pe=group_dada2_filter_trim_pe dada2_dereplicate_fastq_pe=group_dada2_dereplicate_fastq_pe dada2_sample_inference_pe=group_dada2_sample_inference_pe dada2_merge_pairs=group_dada2_merge_pairs dada2_make_table_pe=group_dada2_make_table_pe hisat2_remove_human_sensitive_se=group_hisat2_remove_human_sensitive_se dada2_quality_profile_se=group_dada2_quality_profile_se dada2_filter_trim_se=group_dada2_filter_trim_se dada2_dereplicate_fastq_se=group_dada2_dereplicate_fastq_se dada2_sample_inference_se=group_dada2_sample_inference_se dada2_merge_pairs=group_dada2_merge_pairs dada2_make_table_se=group_dada2_make_table_se get_speciateit_inputs=group_speciateit_valencia speciateit_classify=group_speciateit_valencia speciateit_count_table=group_speciateit_valencia valencia=group_speciateit_valencia adding_csts_to_ps=group_speciateit_valencia get_speciateit_taxa=group_speciateit_valencia \
--group-components group_trim_primers=10 group_hisat2_remove_human_sensitive_pe=10 group_dada2_quality_profile_pe=10 group_dada2_filter_trim_pe=10 group_dada2_dereplicate_fastq_pe=10 group_dada2_sample_inference_pe=10 group_dada2_merge_pairs=10 group_dada2_make_table_pe=10 group_hisat2_remove_human_sensitive_se=10 group_dada2_quality_profile_se=10 group_dada2_filter_trim_se=10 group_dada2_dereplicate_fastq_se=10 group_dada2_sample_inference_se=10 group_dada2_merge_pairs=10 group_dada2_make_table_se=10 group_speciateit_valencia=10;
