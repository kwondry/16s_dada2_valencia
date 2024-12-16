#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=16GB
#SBATCH --partition=short
#SBATCH --job-name=submit_jobs

mkdir -p ./slurm

snakemake --executor slurm --default-resources slurm_partition=short --configfile config.yaml --use-conda --conda-prefix=~/snakemake_conda_prefix --jobname {rulename}.{jobid} --jobs 500 --keep-going --rerun-incomplete