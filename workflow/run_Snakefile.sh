#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G

module purge
module load bioinfo/Snakemake/7.20.0
file="report/dag-$(date +'%d-%m-%Y').pdf"
snakemake --unlock

snakemake --dag | tail -n+3 | dot -Tpdf > $file

#module unload bioinfo/Snakemake/7.20.0

module load devel/Miniconda/Miniconda3
conda config --set channel_priority strict

#module load bioinfo/Snakemake/8.3.1

#snakemake -s Snakefile --workflow-profile /work/project/briefwp3/Adela/GA_CT/config/ -j unlimited --rerun-incomplete --use-conda --use-envmodules --cores 10

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/tools/devel/python/Python-3.11.1/lib/

snakemake -s Snakefile --profile apoublan --rerun-incomplete --use-conda --use-envmodules

# snakemake -s Snakefile --software-deployment-method conda --cores 10 --jobs 100 --rerun-incomplete --use-conda

# can be equally used with cluster modules as: snakemake --use-envmodules (both are included)
# snakemake -s Snakefile --cores 10 --jobs 100 --rerun-incomplete --use-envmodules


# target command line if using profile:
#snakemake --workflow-profile <path> -j unlimited # assuming an unlimited number of jobs \
#--default-resources slurm_account=<account> slurm_partition=<default partition> \
#--configfile config/config.yaml
