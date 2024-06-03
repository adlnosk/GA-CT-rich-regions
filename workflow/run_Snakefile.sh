#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=500M

module purge

module load bioinfo/Snakemake/7.20.0

file="report/dag-$(date +'%d-%m-%Y').pdf"

snakemake --unlock
snakemake --dag | tail -n+3 | dot -Tpdf > $file

module load devel/Miniconda/Miniconda3
conda config --set channel_priority strict
#conda install python=11.1.0


#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/tools/devel/python/Python-3.11.1/lib/


snakemake -s Snakefile --profile ../config --rerun-incomplete --use-conda --use-envmodules --keep-going


# without profile - all is locally
# with newer snake:
# module unload bioinfo/Snakemake/7.20.0
# module load bioinfo/Snakemake/8.3.1
#snakemake -s Snakefile --software-deployment-method conda --cores 10 --jobs 100 --rerun-incomplete --use-conda
