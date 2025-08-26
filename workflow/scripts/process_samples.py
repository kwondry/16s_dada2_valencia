import os
import pandas as pd
import re
import shutil
import argparse
import gzip


def count_lines_in_file(file_path):
    """Count lines in a file, handling both regular and gzipped files."""
    if file_path.endswith('.gz'):
        with gzip.open(file_path, 'rt') as file:
            return sum(1 for line in file)
    else:
        with open(file_path, 'r') as file:
            return sum(1 for line in file)


def is_fastq_file(filename):
    """Check if a file is a FASTQ file (including compressed)."""
    fastq_extensions = ['.fastq', '.fastq.gz', '.fq', '.fq.gz']
    return any(filename.endswith(ext) for ext in fastq_extensions)


def get_sample_id_from_filename(filename, is_paired=False):
    """Extract sample ID from filename, removing FASTQ extensions and paired-end indicators."""
    # Remove common fastq extensions
    sample_id = filename
    for ext in ['.fastq.gz', '.fastq', '.fq.gz', '.fq']:
        if sample_id.endswith(ext):
            sample_id = sample_id[:-len(ext)]
            break
    
    # Remove paired-end indicators if needed
    if is_paired:
        sample_id = sample_id.replace('_R1', '').replace('_R2', '')
    
    return sample_id


def generate_samplesheet(run_sheet, samplesheet_output):

    # outdir
    os.makedirs("outputs/input", exist_ok=True)
    os.makedirs("outputs/logs", exist_ok=True)
    
    #read in csv of runs
    run_sheet = pd.read_csv(run_sheet)

    samplesheet_data = []
    # iterate through csv of runs
    for _, row in run_sheet.iterrows():
        fastq_dir = row['fastq_dir']
        mapping_file = row['mapping_file']
        run_id = row['run_id']
        sequence_type = row['sequence_type']
        region = row['region']

        low_read_samples = []
        for filename in os.listdir(fastq_dir):
            file_path = os.path.join(fastq_dir, filename)
            if os.path.isfile(file_path) and is_fastq_file(filename):
                try:
                    line_count = count_lines_in_file(file_path)
                    
                    #move the file to separate folder if there are no reads so pipeline can run
                    if line_count < 1000:
                        low_read_samples.append(file_path)
                except Exception as e:
                    print(f"Warning: Could not read file {file_path}: {e}")
                    continue

        #creating a .txt log of low read samples for multiqc
        with open("outputs/logs/too_few_reads.txt", 'a') as f:
            for file in low_read_samples:
                f.write(f"{run_id},{file}\n")

        if sequence_type == "single":
            # for each run, extract all samples from the fastq directory and adds it to a sample sheet
            rev_file_path = ''
            for fastq_file in os.listdir(fastq_dir):
                file_path = os.path.join(fastq_dir, fastq_file)
                if is_fastq_file(fastq_file) and not any(file_path in f for f in low_read_samples):
                    sample_id = get_sample_id_from_filename(fastq_file, is_paired=False)
                    fwd_file_path = file_path
                    samplesheet_data.append([sample_id, fwd_file_path, rev_file_path, run_id, mapping_file, sequence_type, region])
                    # print(samplesheet_data)

        if sequence_type == "paired":
            # for each run, extract all samples from the fastq directory and adds it to a sample sheet
            # First collect all R1 and R2 files
            r1_files = {}
            r2_files = {}
            
            for fastq_file in os.listdir(fastq_dir):
                file_path = os.path.join(fastq_dir, fastq_file)
                if is_fastq_file(fastq_file) and not any(file_path in f for f in low_read_samples):
                    if '_R1_' in fastq_file or '_R1.' in fastq_file:
                        sample_id = get_sample_id_from_filename(fastq_file, is_paired=True)
                        r1_files[sample_id] = file_path
                    elif '_R2_' in fastq_file or '_R2.' in fastq_file:
                        sample_id = get_sample_id_from_filename(fastq_file, is_paired=True)
                        r2_files[sample_id] = file_path
            
            # Match R1 and R2 files by sample ID
            for sample_id in r1_files:
                if sample_id in r2_files:
                    fwd_file_path = r1_files[sample_id]
                    rev_file_path = r2_files[sample_id]
                    samplesheet_data.append([sample_id, fwd_file_path, rev_file_path, run_id, mapping_file, sequence_type, region])

    # making samplesheet
    sheet_out = pd.DataFrame(samplesheet_data, columns=['sample_id', 'fwd', 'rev', 'run_id', 'mapping_file', 'sequence_type', 'region'])
    sheet_out.to_csv(samplesheet_output, index=False)


def main(run_sheet, samplesheet_output):

    print("Generating Samplesheet...")
    generate_samplesheet(run_sheet, samplesheet_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("run_sheet")
    parser.add_argument("samplesheet_output")
    args = parser.parse_args()
    
    main(args.run_sheet, args.samplesheet_output)
