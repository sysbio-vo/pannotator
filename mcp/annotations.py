"""Summaries of, and searches in, the GFF3 files Pannotator produces."""

from collections import Counter
from pathlib import Path
from typing import Iterator
from urllib.parse import unquote

from mcp.server.mcpserver.exceptions import ToolError

# `gene` features duplicate their CDS/RNA child, `region` features are the contigs
STRUCTURAL_TYPES = {"gene", "region"}
# characters Bakta percent-encodes in GFF3 attributes
ENCODED_CHARS = set("%;=&,\t\n\r")
SEARCH_FIELDS = {
    "gene": ("gene",),
    "product": ("product", "Name"),
    "locus_tag": ("locus_tag", "ID"),
    "dbxref": ("Dbxref", "Note"),
}


def gff_files(outdir: Path, samples: list[str] | None = None) -> list[Path]:
    if not outdir.is_dir():
        raise ToolError(f"Not an existing directory: {outdir}")
    available = {path.name.removesuffix(".gff3"): path for path in sorted(outdir.glob("*.gff3"))}
    if not available:
        raise ToolError(f"No GFF3 files in {outdir}")
    if samples is None:
        return list(available.values())
    unknown = [sample for sample in samples if sample not in available]
    if unknown:
        raise ToolError(f"No annotation in {outdir} for: {', '.join(unknown[:20])}")
    return [available[sample] for sample in samples]


def gff_lines(path: Path) -> Iterator[str]:
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                return
            yield line.rstrip("\n")


def feature_lines(path: Path) -> Iterator[str]:
    return (line for line in gff_lines(path) if not line.startswith("#") and line.count("\t") == 8)


def parse_attributes(text: str) -> dict[str, str]:
    return {key: unquote(value) for key, _, value in (item.partition("=") for item in text.split(";")) if key}


def is_pseudogene(attributes: dict[str, str]) -> bool:
    return "pseudogene" in attributes or attributes.get("pseudo", "").lower() == "true"


def summarize_gff(path: Path) -> dict:
    types: Counter = Counter()
    contigs = length = hypothetical = pseudogenes = 0
    for line in gff_lines(path):
        # Bakta writes these directives in front of the features of each contig, not at the top of the file
        if line.startswith("##sequence-region"):
            contigs += 1
            length += int(line.split()[3])
        if line.startswith("#") or line.count("\t") != 8:
            continue
        columns = line.split("\t")
        if columns[2] in STRUCTURAL_TYPES:
            continue
        types[columns[2]] += 1
        if columns[2] == "CDS":
            attributes = parse_attributes(columns[8])
            hypothetical += "hypothetical protein" in attributes.get("product", "").lower()
            pseudogenes += is_pseudogene(attributes)
    return {
        "sample": path.name.removesuffix(".gff3"),
        "contigs": contigs,
        "length_bp": length,
        "features": dict(types.most_common()),
        "hypothetical_cds": hypothetical,
        "pseudogenes": pseudogenes,
    }


def summarize(outdir: Path, samples: list[str] | None, limit: int, offset: int) -> dict:
    files = gff_files(outdir, samples)
    return {
        "outdir": str(outdir),
        "n_genomes": len(files),
        "offset": offset,
        "genomes": [summarize_gff(path) for path in files[offset:][:limit]],
    }


def search(
    outdir: Path,
    query: str,
    field: str,
    exact: bool,
    feature_type: str | None,
    samples: list[str] | None,
    limit: int,
) -> dict:
    needle = query.strip().lower()
    if not needle:
        raise ToolError("query is empty")
    # cheap check of the raw line before parsing it; not valid for queries that are encoded in the file
    prefilter = needle if needle.isascii() and not ENCODED_CHARS & set(needle) else ""
    keys = SEARCH_FIELDS.get(field) or tuple(key for group in SEARCH_FIELDS.values() for key in group)

    files = gff_files(outdir, samples)
    hits, hits_per_sample = [], {}
    for path in files:
        sample = path.name.removesuffix(".gff3")
        n_hits = 0
        for line in feature_lines(path):
            if prefilter not in line.lower():
                continue
            seqid, _, kind, start, end, _, strand, _, raw = line.split("\t")
            if kind != feature_type if feature_type else kind in STRUCTURAL_TYPES:
                continue
            attributes = parse_attributes(raw)
            if not any(matches(needle, attributes[key].lower(), key, exact) for key in keys if key in attributes):
                continue
            n_hits += 1
            if len(hits) < limit:
                hit = {"sample": sample, "contig": seqid, "type": kind}
                hit |= {"start": int(start), "end": int(end), "strand": strand}
                for key in ("locus_tag", "gene", "product", "Dbxref", "Note"):
                    if key in attributes:
                        hit[key.lower()] = attributes[key]
                if is_pseudogene(attributes):
                    hit["pseudogene"] = True
                hits.append(hit)
        if n_hits:
            hits_per_sample[sample] = n_hits

    total = sum(hits_per_sample.values())
    return {
        "n_genomes_searched": len(files),
        "n_genomes_with_hits": len(hits_per_sample),
        "n_hits": total,
        "truncated": total > len(hits),
        "genomes_without_hits": [
            path.name.removesuffix(".gff3") for path in files if path.name.removesuffix(".gff3") not in hits_per_sample
        ][:50],
        "hits": hits,
    }


def matches(needle: str, value: str, key: str, exact: bool) -> bool:
    if not exact:
        return needle in value
    if key in ("Dbxref", "Note"):
        # lists such as 'SO:0001217,UniRef:UniRef50_P76484', match with or without the prefix
        return any(needle in (item, item.partition(":")[2]) for item in value.split(","))
    return needle == value
