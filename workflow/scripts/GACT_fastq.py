# GA_CT from fasta file (unzipped)

# sbatch --mem 50G --wrap "module load devel/python/Python-3.7.9; python fa_GA_CT_regions_slide.py"

from Bio import SeqIO
import gzip
input_file_path = snakemake.input.fastq_gzip_file
print(input_file_path)
output_file = snakemake.output.out

window_size = snakemake.params.window_size
threshold = snakemake.params.threshold


enriched_regions = []
# window_size = 50
# threshold = 0.95

with gzip.open(input_file_path, "rt") as inf:
    for n, record in enumerate(SeqIO.parse(inf, "fastq")):
        sequence = record.seq
        i_min = 0
        w = 'will see'
        id = record.name
        for i in range(0, len(sequence) - window_size + 1, window_size // 2):
            start = i_min
            end = i + window_size
            # calculate small window frequency
            window = sequence[i:i+window_size]
            ga_count = window.count('G') + window.count('A')
            ct_count = window.count('C') + window.count('T')
            ga_ratio = ga_count / window_size
            ct_ratio = ct_count / window_size
            if (ga_ratio > threshold or ct_ratio > threshold) and (end < (len(sequence)-window_size)):
                i_min = start
                w = 'yes'
            elif (ga_ratio > threshold or ct_ratio > threshold) and (end > (len(sequence)-window_size)) and (end-start > 10*window_size):
                w = 'yes'
                enriched_regions.append(
                    (id, start, end, end-start, sequence[start:end]))
                # print((id, start, end, end-start, sequence[start:end]))
                i_min = i+(window_size // 2)
            elif w == 'yes' and (ga_ratio < threshold and ct_ratio < threshold) and (end-start > 10*window_size):
                w = 'no'
                enriched_regions.append(
                    (id, start, end, end-start, sequence[start:end]))
                # print((id, start, end, end-start, sequence[start:end]))
                i_min = i+(window_size // 2)
            else:
                w = 'no'
                i_min = i+(window_size // 2)


with open(output_file, "w") as f:
    for a, b, c, d, e in enriched_regions:
        f.write(f"{a}\t{b}\t{c}\t{d}\t{e}\n")
