<div align="center">
<img src="pannotator.svg" alt="Pannotator Logo" width="300"/>
</div>
<br>
<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-darkgreen)](https://github.com/sysbio-vo/pannotator/blob/main/LICENSE.md)
[![Run with Conda](http://img.shields.io/badge/Run%20with-Conda-44A833?logo=anaconda)](https://docs.conda.io/en/latest/)
[![Run with Docker](https://img.shields.io/badge/Run%20with-Docker-2496ED?&logo=docker)](https://www.docker.com/)
[![Run with Singularity](https://img.shields.io/badge/Run%20with-Singularity-1d355c)](https://sylabs.io/docs/)

</div>

# Pannotator: prokaryotic genome annotation _at scale_

Pannotator is a scalable and robust pangenome-based prokaryotic genome annotation tool, designed to efficiently process thousands of genomes. It is built upon [Bakta](https://github.com/oschwengers/bakta) to reliably annotate **protein-coding** and **ncRNA** genes, while leveraging the scalability and reproducibility of [Nextflow](https://www.nextflow.io/).

## Description

- Pannotator orchestrates Bakta annotation steps in a modular Nextflow pipeline. It supports the annotation of ncRNA cis-regulatory regions, oriC/oriV/oriT, assembly gaps, as well as tRNA, tmRNA, rRNA, ncRNA genes, CRISPRs, CDSs, and pseudogenes [via Bakta](https://github.com/oschwengers/bakta?tab=readme-ov-file#description).
- To minimise redundant computation, Pannotator clusters CDS features across genomes and annotates only representative sequences from each cluster, propagating annotations back to cluster members.

## Installation

Prerequisites:

- [Nextflow](https://www.nextflow.io/docs/latest/install.html) `>= 21.04.0`
- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Docker](https://www.docker.com/get-started/)/[Singularity](https://sylabs.io/singularity/)

## Examples

To annotate a folder of genomes using an existing Bakta database:

```bash
nextflow run main.nf --indir /path/to/folder/with/genomes --outdir /path/to/output/folder --bakta_db /path/to/bakta_db
```

To output intermediate files, such as MMseqs2 clustering results, proteome FASTA file containing all unique sequences, and others, use the `--save_intermediate` flag.

```bash
nextflow run main.nf --indir /path/to/folder/with/isolates/ --outdir /path/to/output/folder --bakta_db /path/to/bakta_db --save_intermediate
```

If you don't have a Bakta database, the most recent version will be automatically downloaded during the first run. Note that it might take some time, as the `light` database v6.0 is ~1.3 GB, while the `full` database is ~33.1 GB. The `light` database is downloaded by default. You can specify the required database type through the command line:

```bash
nextflow run main.nf --indir /path/to/folder/with/isolates/ --outdir /path/to/output/folder --bakta_db /path/to/save/bakta/db/ --bakta_db_type [light|full]
```

Available generic execution profiles adapted from the [base config by PaM](https://github.com/sanger-pathogens/nextflow-commons/blob/master/configs/nextflow.config).:

- `standard` (default)
- `docker`
- `singularity`
- `conda`

## Examples for the **Wellcome Sanger Institute's** Farm users

Instead of cloning the repository, you can use a dedicated `pannotator` module on Farm, which is maintained to be up to date with the upstream codebase. To start using it, you need to load the environment first:

```bash
module load PaM/environment
module load pannotator
```

After that, you can use the tool by calling `pannotator`. We recommend using the `sanger_lsf` profile when running the pipeline on Farm. For instance, to annotate a folder of genomes with Pannotator, run:

```bash
pannotator --indir /path/to/folder/with/genomes --outdir /path/to/output/folder --bakta_db /path/to/bakta/db -profile sanger_lsf
```
