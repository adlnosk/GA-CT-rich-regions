# detect GACT regions in all reads

import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import os

# Input and output file paths
input_file_path = snakemake.input.fastq_gzip_file
output_file = snakemake.output.out

# Ensure output directory exists
output_dir = os.path.dirname(output_file)
if output_dir and not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Parameters
window_size = snakemake.params.window_size
threshold = snakemake.params.threshold

# Function to calculate GA/CT ratios
def calculate_ratios(window, window_size):
    ga_count = window.count('G') + window.count('A')
    ct_count = window.count('C') + window.count('T')
    ga_ratio = ga_count / window_size
    ct_ratio = ct_count / window_size
    return ga_ratio, ct_ratio

# Process the FASTQ file with error handling
try:
    with gzip.open(input_file_path, "rt") as inf, open(output_file, "w") as outf:
        for title, sequence, _ in FastqGeneralIterator(inf):  # Skip quality scores
            i_min = 0
            w = "will see"
            contig_id = title.split()[0]  # Use the first part of the title as ID

            for i in range(0, len(sequence) - window_size + 1, window_size // 2):
                start = i_min
                end = i + window_size
                window = sequence[i:i + window_size]

                # Calculate GA/CT ratios
                ga_ratio, ct_ratio = calculate_ratios(window, window_size)

                if (ga_ratio > threshold or ct_ratio > threshold) and (end < (len(sequence) - window_size)):
                    i_min = start
                    w = "yes"
                elif (ga_ratio > threshold or ct_ratio > threshold) and (end > (len(sequence) - window_size)) and (end - start > window_size):
                    w = "yes"
                    outf.write(f"{contig_id}\t{start}\t{end}\t{end - start}\t{sequence[start:end]}\n")
                    i_min = i + (window_size // 2)
                elif w == "yes" and (ga_ratio < threshold and ct_ratio < threshold) and (end - start > window_size):
                    w = "no"
                    outf.write(f"{contig_id}\t{start}\t{end}\t{end - start}\t{sequence[start:end]}\n")
                    i_min = i + (window_size // 2)
                else:
                    w = "no"
                    i_min = i + (window_size // 2)
except EOFError:
    print(f"Error: The input file '{input_file_path}' is corrupted or incomplete.")
    raise

