#!/usr/bin/env python3

import argparse
import hashlib
import re
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from BCBio import GFF
from jinja2 import Environment, FileSystemLoader

LOCATION_REGEX = r"\[(?P<start>\d+):(?P<stop>\d+)\]\((?P<strand>[+-]{1})\)"
QUALIFIER_FIELDS = ["Dbxref", "product", "gene", "source"]


def seq_to_hash(seq: str) -> str:
    return hashlib.sha256(seq.encode("utf-8")).hexdigest()


########################################
# -------- GFF3 PARSING -------------- #
########################################


def parse_gff3(gff_path: str, limit_info: dict | None = None) -> dict:

    features = []

    # limit_info = dict(gff_type=["CDS"]) # example

    with open(gff_path, "r") as in_handle:
        for contig_rec in GFF.parse(in_handle, limit_info=limit_info):
            contig_id = contig_rec.id
            for contig_feature in contig_rec.features:

                if contig_feature.type == "region":
                    continue

                regex_location = re.match(LOCATION_REGEX, str(contig_feature.location))
                if regex_location is None:
                    continue
                start, stop, strand = (
                    regex_location.group("start"),
                    regex_location.group("stop"),
                    regex_location.group("strand"),
                )

                feature_dict = {
                    "seq_id": seq_to_hash(str(contig_feature.location.extract(contig_rec.seq))),
                    "contig_id": contig_id,
                    "type": contig_feature.type,
                    "start": int(start),
                    "stop": int(stop),
                    "strand": strand,
                    "length": int(stop) - int(start),
                }

                for qual_field in QUALIFIER_FIELDS:
                    if qual_field in contig_feature.qualifiers:
                        feature_dict[qual_field] = contig_feature.qualifiers[qual_field]
                    else:
                        feature_dict[qual_field] = None

                features.append(feature_dict)

    return pd.DataFrame(features)


########################################
# -------- PER GENOME STATS ---------- #
########################################


def compute_genome_stats(df):
    stats = {}

    cds = df[df.type == "CDS"]
    rna = df[df.type.str.contains("RNA", na=False)]

    stats["cds_count"] = len(cds)
    stats["rna_count"] = len(rna)
    stats["mean_cds_length"] = cds.length.mean() if len(cds) else 0
    stats["median_cds_length"] = cds.length.median() if len(cds) else 0

    stats["strand_balance"] = df.strand.value_counts(normalize=True).to_dict()

    return stats


########################################
# -------- PANGENOME ----------------- #
########################################


def build_presence_absence(genome_to_df):
    gene_sets = {}

    for genome, df in genome_to_df.items():
        genes = set(df.get("seq_id", []))
        gene_sets[genome] = genes

    all_genes = sorted(set.union(*gene_sets.values()))
    matrix = pd.DataFrame(0, index=genome_to_df.keys(), columns=all_genes)

    for genome, genes in gene_sets.items():
        matrix.loc[genome, list(genes)] = 1

    return matrix


def pangenome_curve(matrix, n_iter=1):
    genomes = list(matrix.index)
    curves = []

    for _ in range(n_iter):
        np.random.shuffle(genomes)
        seen = set()
        counts = []

        for g in genomes:
            genes = set(matrix.columns[matrix.loc[g] == 1])
            seen |= genes
            counts.append(len(seen))

        curves.append(counts)

    return np.mean(curves, axis=0)


########################################
# -------- PLOTTING ------------------ #
########################################


def plot_cds_length(df):
    cds = df[df.type == "CDS"]
    fig = px.histogram(cds, x="length", nbins=50, title="CDS Length Distribution")
    return fig.to_html(full_html=False)


def plot_genome_stats(stats_dict):
    df = pd.DataFrame(stats_dict).T.reset_index().rename(columns={"index": "genome"})
    fig = px.bar(df, x="genome", y="cds_count", title="CDS count per genome")
    return fig.to_html(full_html=False)


def plot_pangenome(curve):
    fig = go.Figure()
    fig.add_trace(go.Scatter(y=curve, mode="lines", name="Pangenome"))
    fig.update_layout(title="Pangenome Saturation Curve", xaxis_title="Genomes", yaxis_title="Genes")
    return fig.to_html(full_html=False)


def plot_presence_absence(matrix):
    fig = px.imshow(matrix, aspect="auto", title="Presence/Absence Heatmap")
    return fig.to_html(full_html=False)


########################################
# MULTIQC-STYLE HTML TEMPLATE
########################################


def render_multiqc_style(outfile, plots, title="Pannotator Report"):
    output_dir = Path(outfile).parent
    output_dir.mkdir(exist_ok=True, parents=True)

    env = Environment(loader=FileSystemLoader("."))

    template = env.from_string(
        """
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{{ title }}</title>

<link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>

<style>
body {
    font-size: 14px;
}

.sidebar {
    height: 100vh;
    position: fixed;
    width: 250px;
    background-color: #1f2933;
    color: white;
    padding: 15px;
}

.sidebar h4 {
    color: #fff;
}

.sidebar a {
    color: #cbd5e1;
    text-decoration: none;
    display: block;
    margin: 8px 0;
}

.sidebar a:hover {
    color: #ffffff;
}

.content {
    margin-left: 260px;
    padding: 20px;
}

.plot-card {
    margin-bottom: 15px;
}

.plot-card .card-body {
    padding: 10px;
}

.plot-container {
    height: 400px;
}

.header-bar {
    background: #111827;
    color: white;
    padding: 10px 20px;
    margin-bottom: 15px;
}

.logo-placeholder {
    width: 40px;
    height: 40px;
    background: #374151;
    display: inline-block;
    margin-right: 10px;
}

</style>
</head>

<body>

<div class="sidebar">
    <h4>Pannotator Report</h4>
    <strong>Pannotator</strong>

    <hr>

    <a href="#overview">Overview</a>
    <a href="#cds">CDS Stats</a>
    <a href="#pangenome">Pangenome</a>
    <a href="#presence">Presence/Absence</a>
</div>

<div class="content">

<div class="header-bar d-flex align-items-center">
    <h5 class="mb-0">{{ title }}</h5>
</div>

<!-- OVERVIEW -->
<section id="overview">
<h4>Overview</h4>
<div class="row">
    {% for name, plot in plots['overview'].items() %}
    <div class="col-lg">
        <div class="card plot-card">
            <div class="card-header">{{ name }}</div>
            <div class="card-body">
                {{ plot | safe }}
            </div>
        </div>
    </div>
    {% endfor %}
</div>
</section>

<!-- CDS -->
<section id="cds" class="mt-4">
<h4>CDS Statistics</h4>
<div class="accordion" id="cdsAccordion">

{% for name, plot in plots['cds'].items() %}
<div class="accordion-item">
    <h2 class="accordion-header">
        <button class="accordion-button collapsed" type="button"
        data-bs-toggle="collapse"
        data-bs-target="#cds_{{ loop.index }}">
            {{ name }}
        </button>
    </h2>
    <div id="cds_{{ loop.index }}" class="accordion-collapse collapse" data-bs-parent="#cdsAccordion">
        <div class="accordion-body">
            {{ plot | safe }}
        </div>
    </div>
</div>
{% endfor %}

</div>
</section>

<!-- PANGENOME -->
<section id="pangenome" class="mt-4">
<h4>Pangenome</h4>
<div class="row">
    {% for name, plot in plots['pangenome'].items() %}
    <div class="col-lg">
        <div class="card plot-card">
            <div class="card-header">{{ name }}</div>
            <div class="card-body">
                {{ plot | safe }}
            </div>
        </div>
    </div>
    {% endfor %}
</div>
</section>

<!-- PRESENCE ABSENCE -->
<section id="presence" class="mt-4">
<h4>Presence / Absence</h4>
<div class="card plot-card">
    <div class="card-header">Heatmap</div>
    <div class="card-body">
        {{ plots['presence']['heatmap'] | safe }}
    </div>
</div>
</section>

</div>
</body>
</html>
"""
    )

    html = template.render(title=title, plots=plots)

    with open(outfile, "w") as f:
        f.write(html)


########################################
# EXPECTED PLOT STRUCTURE
########################################

# plots = {
#   "overview": {
#       "Genome Size vs Genes": html_plot,
#       "CDS Count": html_plot,
#   },
#   "cds": {
#       "CDS Length": html_plot,
#       "Strand Balance": html_plot,
#   },
#   "pangenome": {
#       "Saturation": html_plot,
#       "Gene Frequency": html_plot,
#   },
#   "presence": {
#       "heatmap": html_plot
#   }
# }

########################################
# USAGE
########################################


def main():
    p = argparse.ArgumentParser(description="Generate Pannotator report")
    p.add_argument("--gff3_files", required=True, nargs="+", type=Path, help="GFF3 annotation files")
    p.add_argument("--out", default="pannotator_report.html", type=Path, help="Save path for the report")
    args = p.parse_args()

    gff3_files = args.gff3_files

    genome_to_df = {}
    stats = {}

    for path in gff3_files:
        genome_name = Path(path).stem
        df = parse_gff3(path)
        genome_to_df[genome_name] = df
        stats[genome_name] = compute_genome_stats(df)

    matrix = build_presence_absence(genome_to_df)
    curve = pangenome_curve(matrix)

    plots = {
        "overview": {
            "CDS Count": plot_genome_stats(stats),
        },
        "cds": {"CDS Length": plot_cds_length(list(genome_to_df.values())[0])},
        "pangenome": {
            "Saturation": plot_pangenome(curve),
        },
        "presence": {"heatmap": plot_presence_absence(matrix)},
    }

    render_multiqc_style(args.out, plots)


if __name__ == "__main__":
    main()
