#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G

module purge
module load bioinfo/Snakemake/7.20.0
file="report/dag-$(date +'%d-%m-%Y').pdf"
snakemake --dag | tail -n+2 | dot -Tpdf > $file

module unload bioinfo/Snakemake/7.20.0

module load devel/Miniconda/Miniconda3
conda config --set channel_priority strict

module load bioinfo/Snakemake/8.3.1
snakemake -s Snakefile --software-deployment-method conda --cores 10 --jobs 100 --rerun-incomplete
