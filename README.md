# GA-CT-rich-regions

Workflow to detect assembly regions rich in GA or CT alleles.

### Background

GA/CT simple sequence repeats have been shown to cause a drop in HiFi coverage in the [human genome](https://www.biorxiv.org/content/10.1101/2021.05.26.445798v1) and [tomato genome](https://onlinelibrary.wiley.com/doi/10.1111/tpj.15690). These GA/CT-rich regions can be difficult to assemble, leading to contig breaks. This pipeline aims to detect contigs with problematic regions and gather information to help users link such contig ends together (using multimapped raw reads or HiFi fail reads).

## Workflow structure

- input files: 
	- mandatory: assembly (`.fasta.gz`)
	- optional: HiFi raw reads, HiFi fail reads (`.bam` or `.fastq.gz`)
- main output files:
	- list of GA/CT regions detected in assembly = `file`
```
example table
```
	- list of contigs, that are broken due to GA/CT
```
example table
```
	- list of contig pairs with number of multimapped reads (potential links)
```
example table
```
- additional output files:
	- XX.pdf
	- XX.pdf
	- table for plotting XX
	- table for plotting XX
	

## Repository structure

- `workflow/`: Contains Snakemake pipeline, scripts, environmental modules, and run command.
- `results/`: Directory where the output will be saved.
- `config/`: Contains `config.yaml` to specify cluster environment and `config_path.yaml` to specify input file paths.

## Toy example

DAG of Snakemake rules:
<img src="workflow/report/dag-21-06-2024.pdf" width="500"/>



