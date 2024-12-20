#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=500M

module purge

module load bioinfo/Snakemake/7.20.0

file="report/dag-$(date +'%d-%m-%Y').pdf"

snakemake --unlock
snakemake --dag | tail -n+4 | dot -Tpdf > $file

module load devel/Miniconda/Miniconda3
conda config --set channel_priority strict

snakemake -s Snakefile --profile ../config --rerun-incomplete --use-conda --use-envmodules --keep-going

