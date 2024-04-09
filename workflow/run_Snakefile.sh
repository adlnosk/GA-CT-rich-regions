#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G

module purge
module load bioinfo/Snakemake/7.20.0
snakemake --dag | tail -n+2 | dot -Tpdf > dag.pdf


module purge
module load bioinfo/Snakemake/8.3.1
snakemake -s Snakefile --dryrun
