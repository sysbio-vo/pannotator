<p align="center">
<img src="pannotator.svg" alt="pannotator_logo" width="300"/>
</p>

# Pannotator: prokaryotic genome annotation _at scale_

Pannotator is a scalable and robust pangenome-based prokaryotic genome annotation tool, designed to efficiently process hundreds of genomes. It is built upon [Bakta](https://github.com/oschwengers/bakta) to reliably annotate **protein-coding** and **ncRNA** genes, while leveraging the workflow scalability and reproducibility of [Nextflow](https://www.nextflow.io/).

## Description

- Pannotator orchestrates Bakta annotation steps in a modular Nextflow pipeline. It supports the annotation of ncRNA cis-regulatory regions, oriC/oriV/oriT, assembly gaps, as well as tRNA, tmRNA, rRNA, ncRNA genes, CRISPR, CDS and pseudogenes [via Bakta](https://github.com/oschwengers/bakta?tab=readme-ov-file#description).
- To minimise redundant computation, Pannotator clusters CDS features across genomes and annotates only representative sequences from each cluster, propagating annotations back to cluster members.

## Installation

Prerequisites:

- [Nextflow](https://www.nextflow.io/docs/latest/install.html) `>= 21.04.0`
- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Docker](https://www.docker.com/get-started/)/[Singularity](https://sylabs.io/singularity/)

The Conda environment will be automatically created during the first run. Similarly, if you opt for the `docker` or `singularity` mode, the needed images will be automatically downloaded by Nextflow.

## Examples

Annotate a folder of isolate samples with a full pipeline. This command searches for a Bakta database in the working directory and, if none is found, downloads the `light` Bakta database by default:

```bash
nextflow run main.nf --indir /path/to/folder/with/isolates/
```

Change the output directory with the `--outdir` parameter. If you already have a Bakta database downloaded, pass it as a parameter. For a richer output, save intermediate files with `--save_intermediate`.

A more detailed command may look like this:

```bash
nextflow run main.nf --indir /path/to/folder/with/isolates/ --outdir test_demo_run --save_intermediate --bakta_db /path/to/Bakta/db/ -resume
```
