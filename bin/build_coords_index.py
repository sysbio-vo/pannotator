#!/usr/bin/env python3

from pathlib import Path
from typing import Iterable
import hashlib
import argparse
import json

from Bio import SeqIO
from BCBio import GFF


def find_paired_files(input_dir: Path) -> tuple:
    gffs = {f.name.split('.')[0]: f for f in input_dir.glob(f"*.gff3")}
    faas = {f.name.split('.')[0]: f for f in input_dir.glob(f"*.faa")}
    pairs = {}
    for sample, gff in gffs.items():
        faa = faas.get(sample)
        if faa:
            pairs[sample] = (gff, faa)
    return pairs

def parse_faa(faa_path: Path) -> dict:
    result = {}
    with open(faa_path, 'rt') as f:
        for rec in SeqIO.parse(f, 'fasta'):
            result[rec.id] = str(rec.seq)
    return result

def parse_gff(gff_path: Path) -> Iterable[dict]:
    with open(gff_path, 'rt') as f:
        for rec in GFF.parse(f):
            contig = rec.id
            for f in rec.features:
                if f.type != 'CDS':
                    continue
                locus_tag = f.qualifiers.get('locus_tag')[0]
                coords = f.location
                start = int(coords.start)
                end = int(coords.end)
                strand = '+' if f.location.strand == 1 else '-'
                yield {
                    'contig': contig,
                    'locus_tag': locus_tag,
                    'start': start,
                    'end': end,
                    'strand': strand}
                
def build_index(input_dir: Path) -> dict:
    index = {}
    pairs = find_paired_files(input_dir)
    for sample, (gff, faa) in pairs.items():
        faa_dict = parse_faa(faa)
        for cds in parse_gff(gff):
            locus_tag = cds['locus_tag']
            aa_seq = faa_dict.get(locus_tag)
            if aa_seq is None:
                continue
            aa_seq_hash = hashlib.md5(aa_seq.encode('utf-8')).hexdigest()
            new_member = index.setdefault(aa_seq_hash, {'aa_seq': aa_seq, 'members': []})  
            new_member['members'].append({
                'sample': sample,
                'contig': cds['contig'],
                'locus': locus_tag,
                'start': cds['start'],
                'end': cds['end'],
                'strand': cds['strand'],})
    return index

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build CDS-coords index from cds-only bakta results')
    parser.add_argument('input_dir', type=Path, help='Folder with *.cds-only.gff3 and *.cds-only.faa files')
    parser.add_argument('output_json', type=Path, help='Path to write resulting JSON file')
    args = parser.parse_args()
    print(f"Building CDS-coords index from files in {args.input_dir}")
    print(f"Writing results to {args.output_json}")
    
    index = build_index(args.input_dir)
    args.output_json.write_text(json.dumps(index, indent=2))
    