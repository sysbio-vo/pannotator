#!/usr/bin/env python3

import json
from pathlib import Path
import argparse


def load_annotations(annotation_path):
    with open(annotation_path, 'r') as f:
        return json.load(f)

def load_cds_index(cds_index_path):
    with open(cds_index_path, 'r') as f:
        data = json.load(f)
    cds_index = {}
    for hexdigest, features in data.items():
        for member in features['members']:
            cds_index[(member['sample'], member['locus'])] = hexdigest
    return cds_index


def annotate_gff3(gff3_path, sample, annotations, cds_index, out_path):
    with open (gff3_path, 'r') as f_in, open (out_path, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                f_out.write(line)
                continue
        
            parts = line.strip().split('\t')
            if len(parts) != 9 or parts[2] != "CDS":
                f_out.write(line) 
                continue

            attr = parts[-1]

            attr_dict = dict(item.split('=', 1) for item in attr.split(';') if '=' in item)

            locus = attr_dict.get("locus_tag") 
            hexdigest = cds_index.get((sample, locus))

            if hexdigest:
                annotation = annotations.get(hexdigest)
                if annotation:
                    attr_dict["product"] = annotation.get("product", "hypothetical protein")
                    attr_dict["Name"] = annotation.get("product", "hypothetical protein")

                    if annotation.get("gene"):
                        attr_dict["gene"] = annotation["gene"]

                    if annotation.get("db_xrefs"):
                        attr_dict["Dbxref"] = ",".join(annotation["db_xrefs"])
            else:
                # how will we mark unannotated features?
                attr_dict["unannotated"] = "true"

            new_attr = ';'.join(f"{key}={value}" for key, value in attr_dict.items())
            parts[-1] = new_attr

            f_out.write('\t'.join(parts) + '\n')


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Map CDS annotations with CDS coordinates and produce annotated GFF3 files")
    p.add_argument("--gff3_dir", required=True, help="Directory with *.gff3 files")
    p.add_argument("--cds_index", required=True, help="cds_index.json file")
    p.add_argument("--annotations", required=True, help="bulk_protein_annotations.json file")
    p.add_argument("--out", default="annotated_gff3", help="Output directory")
    args = p.parse_args()

    gff3_dir = Path(args.gff3_dir)
    out_dir = Path(args.out)

    out_dir.mkdir(exist_ok=True, parents=True)

    annotations = load_annotations(args.annotations)
    cds_index = load_cds_index(args.cds_index)

    for gff3_path in gff3_dir.glob("*.gff3"):
        sample = gff3_path.stem.split(".")[0]
        out_path = out_dir / f"{sample}.gff3"
        annotate_gff3(gff3_path, sample, annotations, cds_index, out_path)