<div align="center">
<img src="pannotator.svg" alt="Pannotator Logo" width="300"/>
</div>

<div style="height: 20px;"></div>

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-darkgreen)](https://github.com/sysbio-vo/pannotator/blob/main/LICENSE.md)
[![Run with Conda](http://img.shields.io/badge/Run%20with-Conda-44A833?logo=anaconda)](https://docs.conda.io/en/latest/)
[![Run with Docker](https://img.shields.io/badge/Run%20with-Docker-2496ED?&logo=docker)](https://www.docker.com/)
[![Run with Singularity](https://img.shields.io/badge/Run%20with-Singularity-1d355c)](https://sylabs.io/docs/)

</div>

# Pannotator: prokaryotic genome annotation _at scale_

Pannotator is a scalable and robust pangenome-based prokaryotic genome annotation tool, designed to efficiently process hundreds of genomes. It is built upon [Bakta](https://github.com/oschwengers/bakta) to reliably annotate **protein-coding** and **ncRNA** genes, while leveraging the workflow scalability and reproducibility of [Nextflow](https://www.nextflow.io/).

## Description

- Pannotator orchestrates Bakta annotation steps in a modular Nextflow pipeline. It supports the annotation of ncRNA cis-regulatory regions, oriC/oriV/oriT, assembly gaps, as well as tRNA, tmRNA, rRNA, ncRNA genes, CRISPR, CDS and pseudogenes [via Bakta](https://github.com/oschwengers/bakta?tab=readme-ov-file#description).
- To minimise redundant computation, Pannotator clusters CDS features across genomes and annotates only representative sequences from each cluster, propagating annotations back to cluster members.

## Installation

Prerequisites:

- [Nextflow](https://www.nextflow.io/docs/latest/install.html) `>= 21.04.0`
- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Docker](https://www.docker.com/get-started/)/[Singularity](https://sylabs.io/singularity/)

## Examples

Annotate a folder of isolate samples with a full pipeline. During the first run, the pipeline searches for a Bakta database in the working directory and, if none is found, downloads the `light` Bakta database by default:

```bash
nextflow run main.nf --indir /path/to/folder/with/isolates/ -profile local
```

Change the output directory with the `--outdir` parameter.

```bash
nextflow run main.nf --indir /path/to/folder/with/isolates/ -profile local --outdir test_run
```

For a richer output, save intermediate files with `--save_intermediate`

```bash
nextflow run main.nf --indir /path/to/folder/with/isolates/ -profile local --outdir test_run --save_intermediate
```

If you already have a Bakta database downloaded, pass it as a parameter. By default, the database is assumed to be of type `light`. Make sure to indicate the correct type if needed. This is required to run the annotation steps that rely on the full database, such as pseudogene search.

```bash
nextflow run main.nf --indir /path/to/folder/with/isolates/ -profile local --outdir test_run --save_intermediate --bakta_db /path/to/full/Bakta/db/ --bakta_db_type full
```

Select among other execution profiles.

- `standard` (default)
- `docker`
- `singularity`
- `conda`

For more information regarding the profiles, please refer to the [base config by PaM](https://github.com/sanger-pathogens/nextflow-commons/blob/master/configs/nextflow.config).
