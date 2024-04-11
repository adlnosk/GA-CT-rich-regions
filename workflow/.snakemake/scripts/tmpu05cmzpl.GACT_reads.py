
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/usr/local/bioinfo/src/Snakemake/snakemake-v8.3.1/snakemake-v8.3.1_venv/lib/python3.11/site-packages', '/work/project/briefwp3/Adela/GA_CT/workflow', '/usr/local/bioinfo/src/Snakemake/snakemake-v8.3.1/snakemake-v8.3.1_venv/bin', '/tools/devel/python/Python-3.11.1/lib/python3.11', '/tools/devel/python/Python-3.11.1/lib/python3.11/lib-dynload', '/usr/local/bioinfo/src/Snakemake/snakemake-v8.3.1/snakemake-v8.3.1_venv/lib/python3.11/site-packages', '/home/apoublan/.cache/snakemake/snakemake/source-cache/runtime-cache/tmpt7i49nwp/file/work/project/briefwp3/Adela/GA_CT/workflow/scripts', '/work/project/briefwp3/Adela/GA_CT/workflow/scripts']); import pickle; snakemake = pickle.loads(b'\x80\x04\x954\x0c\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94\x8cQ/work/project/briefwp3/Adela/Medicago_sativa/fastq_file/Medicago_sativa1.fastq.gz\x94a}\x94(\x8c\x06_names\x94}\x94\x8c\x0ffastq_gzip_file\x94K\x00N\x86\x94s\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x12\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x18)}\x94\x8c\x05_name\x94h\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bh\x0eh\nub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94\x8cM/work/project/briefwp3/Adela/GA_CT/results/Medicago_sativa_GACT_raw_reads.txt\x94a}\x94(h\x0c}\x94\x8c\x03out\x94K\x00N\x86\x94sh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bh)h&ub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94(K2G?\xeeffffffe}\x94(h\x0c}\x94(\x8c\x0bwindow_size\x94K\x00N\x86\x94\x8c\tthreshold\x94K\x01N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bh:K2h<G?\xeeffffffub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94\x8c\x0fMedicago_sativa\x94a}\x94(h\x0c}\x94\x8c\x04spec\x94K\x00N\x86\x94sh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94b\x8c\x04spec\x94hKub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01\x8c\x04/tmp\x94e}\x94(h\x0c}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bhbK\x01hdK\x01hfh_ub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bub\x8c\x06config\x94}\x94(\x8c\x07species\x94]\x94\x8c\x0fMedicago_sativa\x94a\x8c\traw_reads\x94}\x94(\x8c\x0fMedicago_sativa\x94\x8c\x97/save/project/pepr_agroeco/WP1-WP2/FORAGE/AGRODIV_MESAT_MILKYMAX/r84087_20231102_141650/1_A01/hifi_reads/m84087_231102_142501_s1.hifi_reads.default.bam\x94\x8c\x14Diplocyclos_palmidus\x94\x8c\x93/save/project/pepr_agroeco/WP1-WP2/FORAGE/AGRODIV_DIPAL_****/r84087_20231011_153509/1_C01/hifi_reads/m84087_231011_164406_s3.hifi_reads.default.bam\x94\x8c\x12Dactylis_glomerata\x94\x8c\x8c/save/project/pepr_agroeco/WP1-WP2/FORAGE/AGRODIV_DAGLO_LUKIR/r84087_20231220_162434/1_C01/hifi_reads/m84087_231220_203121_s3.hifi_reads.bam\x94u\x8c\x08assembly\x94}\x94(\x8c\x0fMedicago_sativa\x94\x8cU/work/project/briefwp3/Adela/Medicago_sativa/assembly/Medicago_sativa.asm.bp.p_ctg.fa\x94\x8c\x14Diplocyclos_palmidus\x94\x8cf/work/project/briefwp3/Adela/Diplocyclos_palmidus/downsampled_ass/Diplocyclos_palmidus.asm.bp.p_ctg.fa\x94\x8c\x12Dactylis_glomerata\x94\x8c[/work/project/briefwp3/Adela/Dactylis_glomerata/assembly/Dactylis_glomerata.asm.bp.p_ctg.fa\x94u\x8c\x05fastq\x94}\x94(\x8c\x0fMedicago_sativa\x94\x8cQ/work/project/briefwp3/Adela/Medicago_sativa/fastq_file/Medicago_sativa1.fastq.gz\x94\x8c\x14Diplocyclos_palmidus\x94\x8cZ/work/project/briefwp3/Adela/Diplocyclos_palmidus/fastq_file/Diplocyclos_palmidus.fastq.gz\x94\x8c\x12Dactylis_glomerata\x94\x8cV/work/project/briefwp3/Adela/Dactylis_glomerata/fastq_file/Dactylis_glomerata.fastq.gz\x94u\x8c\nfail_reads\x94}\x94(\x8c\x0fMedicago_sativa\x94\x8c\x97/save/project/pepr_agroeco/WP1-WP2/FORAGE/AGRODIV_MESAT_MILKYMAX/r84087_20231102_141650/1_A01/fail_reads/m84087_231102_142501_s1.fail_reads.default.bam\x94\x8c\x14Diplocyclos_palmidus\x94\x8c\x95/save/project/pepr_agroeco/WP1-WP2/FORAGE/AGRODIV_DIPAL_GENOME/r84087_20231011_153509/1_C01/fail_reads/m84087_231011_164406_s3.fail_reads.default.bam\x94\x8c\x12Dactylis_glomerata\x94\x8c\x8c/save/project/pepr_agroeco/WP1-WP2/FORAGE/AGRODIV_DAGLO_LUKIR/r84087_20231220_162434/1_C01/fail_reads/m84087_231220_203121_s3.fail_reads.bam\x94u\x8c\rGACT_assembly\x94\x8c\x04True\x94\x8c\tGACT_fail\x94\x8c\x04True\x94\x8c\x10GACT_multimapped\x94\x8c\x04True\x94\x8c\rGACT_bridging\x94\x8c\x04True\x94u\x8c\x04rule\x94\x8c\x0bGA_CT_reads\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c3/work/project/briefwp3/Adela/GA_CT/workflow/scripts\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/work/project/briefwp3/Adela/GA_CT/workflow/scripts/GACT_reads.py';
######## snakemake preamble end #########
# GA_CT from fasta file (unzipped)

# sbatch --mem 50G --wrap "module load devel/python/Python-3.7.9; python fa_GA_CT_regions_slide.py"

import sys
print(sys.version)

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
