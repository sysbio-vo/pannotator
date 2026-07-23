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

# Pannotator

Pannotator is a scalable pangenome-based tool for prokaryotic genome annotation, built to process collections of thousands of genomes. It wraps [Bakta](https://github.com/oschwengers/bakta) for annotation of protein-coding and ncRNA genes within a [Nextflow](https://www.nextflow.io/) workflow, providing portability and reproducible execution across compute environments. Rather than annotating each genome independently, Pannotator predicts coding sequences across the whole collection, clusters them by sequence identity, annotates a single representative per cluster, and propagates the result to all members, reducing redundant computation from per-genome to per-cluster while keeping annotations consistent across the pangenome.

Pannotator is presented in more detail in our accompanying paper (TBD).

## Installation

Prerequisites:

- [Nextflow](https://www.nextflow.io/docs/latest/install.html) `>= 21.04.0`
- One of [Conda](https://docs.conda.io/), [Docker](https://www.docker.com/), or [Singularity](https://sylabs.io/singularity/) to supply the Bakta and MMseqs2 tools

The pipeline pulls container images (Bakta and MMseqs2) automatically for the container-based profiles. The `conda` profile builds an environment from `environment/pannotator.yaml`. No manual tool installation is required.

Clone the repository with its submodules and add the `nxf_bin` directory to your `PATH`:

```bash
git clone --recurse-submodules https://github.com/sysbio-vo/pannotator.git
cd pannotator
echo "export PATH=\"\$PATH:$(pwd)/nxf_bin\"" >> ~/.bashrc
source ~/.bashrc
```

This makes the `pannotator` command available. Run `pannotator --help` to list parameters.

## Quick start

If no Bakta database is found at `--bakta_db`, a database of the type given by `--bakta_db_type` (`light` by default) is downloaded automatically on the first run. This can take a while: the `light` database is ~1.3 GB, the `full` database ~33.1 GB.

```bash
pannotator --indir /path/to/assemblies/ --outdir /path/to/output/ -profile standard
```

## Input requirements

- `--indir` must be a directory containing assemblies in FASTA format (optionally gzipped).
- By default every file in the directory is treated as an input. Restrict this with `--infile_extension` (for example `--infile_extension .fasta`).

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--indir` | (required) | Directory of input FASTA assemblies. |
| `--infile_extension` | `""` | Only files ending in this string are used as input. Empty means all files. |
| `--outdir` | `./pannotator_results` | Output directory. |
| `--save_intermediate` | `false` | Publish intermediate files (per-genome CDS/RNA predictions, clustering files, raw annotations) to `--outdir`. |
| `--scratch` | `false` | Run Bakta processes in a node-local scratch directory. Useful for some (cloud) executors. |
| `--bakta_db` | `./bakta_db/db-light` or `./bakta_db/db` | Path to an existing Bakta database. If it does not exist, the database is downloaded here. |
| `--bakta_db_type` | `light` | Bakta database type, `light` or `full`. Pseudogene detection runs only with `full`. |
| `--bakta_args` | `""` | Extra arguments passed through to Bakta. |
| `--compliant` | `false` | Produce INSDC-compliant output (passes `--compliant` to Bakta). |
| `--user_proteins` | `null` | FASTA of expert proteins for CDS annotation (Bakta `--proteins`). |
| `--user_hmms` | `null` | HMMER profile file for CDS annotation (Bakta `--hmms`). |
| `--auxiliary_db` | `null` | Path to a pangenome annotation index (JSON). See [Reusing annotations](#reusing-annotations-across-runs). |
| `--extend_auxdb` | `false` | Add annotations from the current run to an existing pangenome index. |
| `--mmseqs_command` | `easy-linclust` | MMseqs2 clustering command, `easy-linclust` or `easy-cluster`. |
| `--mmseqs_args` | `--min-seq-id 1.0 -c 1.0 --alignment-mode 3` | Arguments passed to the MMseqs2 clustering command (see [Clustering](#clustering)). |
| `--help` | | Print help and exit. |

Nextflow `-profile` selects the execution and container environment. The default `standard` profile runs locally. Container and environment profiles `docker`, `singularity`, and `conda` are adapted from the [base config by PaM](https://github.com/sanger-pathogens/nextflow-commons/blob/master/configs/nextflow.config). LSF profiles (`sanger_lsf`, `conda_lsf`) are also provided.

## Clustering

Coding sequences from all genomes are pooled into a single FASTA file and clustered with MMseqs2 before annotation. Clustering behaviour is controlled by two parameters: `--mmseqs_command` selects the algorithm, and `--mmseqs_args` is passed straight to that command.

Clustering stringency is governed by two MMseqs2 options:

- `--min-seq-id` sets the minimum pairwise sequence identity for two sequences to be placed in the same cluster. Use `0.90` for 90% identity and `0.95` for 95% identity. The default is `1.0`, so only identical sequences cluster together.
- `-c` sets the required alignment coverage. The default here is `1.0` (full-length match), which pairs with the default 100%-identity clustering. When clustering below 100% identity, a coverage of `0.8` is the common choice, following UniRef-style clustering.

When clustering below 100% identity, the representative's annotation is propagated to every cluster member, so members that are similar but not identical to the representative still receive its annotation.

`easy-linclust` runs linear-time clustering: it is faster and uses less memory, which suits large input sets, but it is less sensitive and may miss some distant matches. `easy-cluster` is more sensitive because it compares more sequence pairs, at higher runtime and memory cost. The default here is `easy-linclust`. For large collections keep `easy-linclust`; switch to `easy-cluster` when sensitivity matters more than speed.

90% identity clustering:

```bash
pannotator \
    --indir /path/to/assemblies/ \
    --outdir /path/to/output/ \
    -profile standard \
    --mmseqs_command easy-linclust \
    --mmseqs_args "--min-seq-id 0.9 -c 0.8 --alignment-mode 3"
```

95% identity clustering:

```bash
pannotator \
    --indir /path/to/assemblies/ \
    --outdir /path/to/output/ \
    -profile standard \
    --mmseqs_command easy-linclust \
    --mmseqs_args "--min-seq-id 0.95 -c 0.8 --alignment-mode 3"
```

## Reusing annotations across runs

Annotation of representative proteins is the most expensive step. To avoid repeating it, save the annotations to a pangenome index and reuse them in later runs with `--auxiliary_db`:

```bash
pannotator --indir /path/to/assemblies/ --outdir /path/to/output/ \
    --bakta_db /path/to/bakta_db --auxiliary_db /path/to/pangenome_index.json
```

If the index already exists, CDS annotations are fetched from it and only proteins absent from the index are annotated with Bakta. To add the current run's annotations to an existing index, add `--extend_auxdb`:

```bash
pannotator --indir /path/to/assemblies/ --outdir /path/to/output/ \
    --bakta_db /path/to/bakta_db --auxiliary_db /path/to/pangenome_index.json --extend_auxdb
```

## Output files

Outputs are written under `--outdir` (default `./pannotator_results`). The per-genome GFF3 files are the main result and are always published. The other files are intermediates, published only when `--save_intermediate` is `true`.

| Path | Description |
|------|-------------|
| `<sample>.gff3` | Final per-genome annotation (CDS, RNA features, and short-ORF / extra search combined). Always published. |
| `annotated_pkl/*.pkl` | Per-genome CDS features with cluster annotations merged in. |
| `CDSS_bakta/<sample>.cds-only.faa`, `.cds-only.pkl` | Per-genome CDS prediction. |
| `RNAS_bakta/<sample>.rna-only.pkl` | Per-genome RNA prediction. |
| `mmseqs_clustering_all_seqs.fasta`, `mmseqs_clustering_cluster.tsv`, `mmseqs_clustering_rep_seq.fasta` | MMseqs2 clustering: all sequences, cluster membership, and cluster representatives. |
| `annotated_proteins_bakta/unique_proteins_annotation.json` | Bakta annotation of the representative proteins. |
| `bulk_protein_annotations.json` | Representative-protein annotations reduced to a lookup keyed by protein hash. |
| `bulk_protein_annotations_extended.json` | Annotation lookup extended to cluster members (sub-100% identity clustering only). |
| `cds_with_pseudogenes/*with_pseudogenes.pkl` | CDS features with pseudogenes, produced only when `--bakta_db_type full`. |

When `--auxiliary_db` is used, the pangenome index JSON is written next to the path given by `--auxiliary_db`.

## How the pipeline works

1. Bakta database setup. If `--bakta_db` does not exist, the database of type `--bakta_db_type` is downloaded.
2. CDS prediction. Bakta runs per genome in `--cds-only` mode, producing protein FASTA and a pickled feature object for each assembly.
3. Clustering. Proteins from all genomes are pooled and clustered with MMseqs2, yielding one representative sequence per cluster.
4. Representative annotation. Bakta annotates only the representative proteins (`bakta_proteins`). With `--auxiliary_db`, known proteins are taken from the index and only the remaining ones are annotated. The result is reduced to a per-protein lookup.
5. Annotation propagation and merge. For sub-100% identity clustering, annotations are propagated from representatives to all cluster members. Each genome's CDS features are then updated with the annotation of their cluster representative.
6. Pseudogene detection. With `--bakta_db_type full`, pseudogenes are detected across the annotated CDS set (`bakta_pseudo_bulk`). Skipped for the `light` database.
7. RNA prediction. Bakta runs per genome in `--rna-only` mode.
8. Short-ORF / extra search. For each genome, Bakta combines the assembly, its annotated CDS, and its RNA features in `--sorf-extra` mode to produce the final per-genome GFF3.
9. Pangenome index update. With `--auxiliary_db` (new index or `--extend_auxdb`), the run's annotations are written back to the index for reuse.

## Citation

If you use Pannotator, please cite both Pannotator and Bakta. Citation metadata for Pannotator is in [CITATION.cff](CITATION.cff).
The accompanying paper is TBD.

## Credits and Affiliations

This tool is developed and maintained in collaboration between the following research institutions:

- **Institute of Molecular Biology and Genetics of the NAS of Ukraine** – [Systems Biology Group](https://sysbio.org.ua/team)
- **Wellcome Sanger Institute** – [Parasites and Microbes Informatics](https://www.sanger.ac.uk/group/parasites-and-microbes-informatics/)

**Funding:** This work was supported by the Wellcome Trust.

## License

MIT. Copyright © 2024 Genome Research Ltd. See [LICENSE.md](LICENSE.md).
