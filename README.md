# 16s_dada2_valencia

A Snakemake workflow for processing 16S rRNA sequencing data using DADA2 and VALENCIA for community state type (CST) classification.

## Documentation

For detailed documentation, visit: https://kwondry.github.io/documentation/materials/16S/valencia/

## Prerequisites

1. **Required Software**:
   - Snakemake (v8.20 or higher)
   - Conda or Mamba
   - Python 3.10 or higher
   - R 4.3 or higher

2. **Required Databases**:
   - HISAT2 human reference database
   - GTDB reference databases
   - SpeciateIT VALENCIA database

## Installation

1. Install Snakemake and Snakedeploy:
   ```bash
   conda install -c conda-forge -c bioconda snakemake snakedeploy
   ```

2. Create a new project directory and deploy the workflow:
   ```bash
   mkdir -p 16s_dada2_valencia
   cd 16s_dada2_valencia
   snakedeploy deploy-workflow https://github.com/kwondry/16s_dada2_valencia . --branch main
   ```

## Configuration

1. **Database Paths**:
   Edit `config/config.yaml` to set the paths to your reference databases:
   ```yaml
   hisat_db: "/path/to/hisat2/human/reference"  # Path to HISAT2 index of T2T human genome
   gtdb_tax_db: "/path/to/GTDB/taxonomy/database"  # DADA2 GTDB taxonomy database
   gtdb_species_db: "/path/to/GTDB/species/database"  # DADA2 GTDB species database 
   speciateit_db: "/path/to/speciateit_valencia"  # SpeciateIT VALENCIA database
   ```

   To obtain these databases:
   - **HISAT2 human reference**: Download the T2T human genome from NCBI (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/) and index it using HISAT2's `hisat2-build` command
   - **GTDB databases**: Download the DADA2-formatted GTDB reference files from Zenodo (https://zenodo.org/record/4587955)
   - **SpeciateIT database**: Clone the VALENCIA database from the Ravel Lab GitHub (https://github.com/ravel-lab/speciateIT)

2. **DADA2 Parameters**:
   The workflow includes default parameters for DADA2 processing. You can modify these in `config/config.yaml`:
   ```yaml
   dada2:
     pe:
       maxEE: 3
       truncLen: [286, 260]
       trimLeft: [10, 10]
       merge:
         minOverlap: 8
         maxMismatch: 1
     se:
       maxEE: 1
       truncLen: 230
       trimLeft: 10
   ```

3. **Primer Sequences**:
   Update the primer sequences in `config/config.yaml` if using different primers:
   ```yaml
   primers:
     v3v4:
       forward: "ACTCCTRCGGGAGGCAGCAG"
       reverse: "GGACTACHVGGGTWTCTAAT"
   ```

## Input Data

1. Create a sequencing runsheet CSV file with the following columns:
   - `run_id`: Run identifier
   - `mapping_file`: Path to mapping file containing sample IDs
   - `fastq_dir`: Path to directory containing FASTQ files
   - `region`: Sequencing region used (v4 or v3v4)

2. The mapping file should contain sample identifiers that match the FASTQ filenames in the fastq_dir.

3. The workflow will automatically process both single-end and paired-end data based on the FASTQ files present.

## Running the Workflow

1. **Test Run**:
   ```bash
   snakemake -n
   ```

2. **Local Execution**:
   ```bash
   snakemake --use-conda --cores all
   ```

3. **SLURM Execution**:
   ```bash
   snakemake \
       --executor slurm \
       --default-resources slurm_partition=short \
       --jobs 100 \
       --use-conda \
       --rerun-incomplete
   ```

## Output

The workflow generates two main outputs in the `outputs/results` directory:
- A phyloseq object containing ASV counts, taxonomy assignments, CST assignment, and sample metadata
- A MultiQC report summarizing read quality metrics and taxonomic composition


![image](https://github.com/user-attachments/assets/4a850cce-ee16-4d00-97e1-f1616e9bb0dd)

