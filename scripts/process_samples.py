import os
import pandas as pd
import re
import shutil
import argparse


def generate_samplesheet(run_sheet, samplesheet_output):

    # outdir
    os.makedirs("outputs/input", exist_ok=True)

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

        # moves samples with 0 reads into a separate folder
        print("Checking if samples have 0 reads")
        no_reads_dir = os.path.join(fastq_dir, "no_reads")
        if not os.path.exists(no_reads_dir):
            os.makedirs(no_reads_dir)

        for filename in os.listdir(fastq_dir):
            file_path = os.path.join(fastq_dir, filename)
            if os.path.isfile(file_path):
                with open(file_path, 'r') as file:
                    line_count = sum(1 for line in file)
                
                #move the file to separate folder if there are no reads so pipelin can run
                if line_count == 0:
                    if not os.path.exists(no_reads_dir):
                        os.makedirs(no_reads_dir)
                    shutil.move(file_path, no_reads_dir)
                    print(f"Moved {file_path} to '{no_reads_dir}'")

        if sequence_type == "single":
            # for each run, extract all samples from the fastq directory and adds it to a sample sheet
            rev_file_path = ''
            for fastq_file in os.listdir(fastq_dir):
                if fastq_file.endswith('.fastq'):
                    sample_id = fastq_file.replace('.fastq', "")
                    fwd_file_path = os.path.join(fastq_dir, fastq_file)
                    samplesheet_data.append([sample_id, fwd_file_path, rev_file_path, run_id, mapping_file, sequence_type, region])
                    # print(samplesheet_data)

        if sequence_type == "paired":
            # for each run, extract all samples from the fastq directory and adds it to a sample sheet
            for fastq_file in os.listdir(fastq_dir):
                if fastq_file.endswith('.fastq') and ('_R1_' in fastq_file or '_R1.' in fastq_file):
                    fwd_file_path = os.path.join(fastq_dir, fastq_file)
                if fastq_file.endswith('.fastq') and ('_R2_' in fastq_file or '_R2.' in fastq_file):
                    sample_id = fastq_file.replace('_R1', "").replace('_R2', "").replace('.fastq', "")
                    rev_file_path = os.path.join(fastq_dir, fastq_file)
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