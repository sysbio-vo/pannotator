import argparse
from pathlib import Path

import pandas as pd
from BCBio import GFF
from collections import defaultdict


def parse_gff3_to_dict(gff_path):
    features = {}

    with open(gff_path, 'r') as gff_file:
        for record in GFF.parse(gff_file):
            contig_id = record.id

            for feature in record.features:
                location = str(feature.location)
                feature_key = f"{contig_id}:{location}"

                qualifiers = {}
                for qual_key, qual_value in feature.qualifiers.items():
                    if qual_key == 'ID' or qual_key == 'locus_tag':
                        continue

                    if isinstance(qual_value, list):
                        qualifiers[qual_key] = ','.join(map(str, qual_value))
                    else:
                        qualifiers[qual_key] = str(qual_value)
                
                features[feature_key] = {
                    'contig_id': contig_id,
                    'type': feature.type,
                    'location': location,
                    **qualifiers
                }

    return features


def compare_gff3_files(gff_ref_path, gff_test_path):
    features_ref = parse_gff3_to_dict(gff_ref_path)
    features_test = parse_gff3_to_dict(gff_test_path)

    keys_ref = features_ref.keys()
    keys_test = features_test.keys()
    
    print(f"Total features in {gff_ref_path.name}: {len(keys_ref)}")
    print(f"Total features in {gff_test_path.name}: {len(keys_test)}")
    print(f"Common features: {len(keys_ref & keys_test)}")
    print(f"Only in {gff_ref_path.name}: {len(keys_ref - keys_test)}")
    print(f"Only in {gff_test_path.name}: {len(keys_test - keys_ref)}")
    print()
    print('-'*80)

    # features only in ref file
    if keys_ref - keys_test:
        print(f"Features only in {gff_ref_path.name}:\n")
        for key in sorted(keys_ref - keys_test):
            f = features_ref[key]
            print(f"Location: {key}")
            print(f"  Type: {f['type']}")
            if 'ID' in f:
                print(f"  ID: {f['ID']}")
            if 'Name' in f:
                print(f"  Name: {f['Name']}")
            print()
    
    print('-'*80)

    # features only in test file
    if keys_test - keys_ref:
        print(f"Features only in {gff_test_path.name}:")
        for key in sorted(keys_test - keys_ref):
            f = features_test[key]
            print(f"Location: {key}")
            print(f"  Type: {f['type']}")
            if 'ID' in f:
                print(f"  ID: {f['ID']}")
            if 'Name' in f:
                print(f"  Name: {f['Name']}")
            print()
    
    print('-'*80)

    # compare common features
    common_keys = keys_ref & keys_test
    if common_keys:
        print(f"Difference in common features:\n")
        differences = defaultdict(list)

        for key in common_keys:
            f_ref = features_ref[key]
            f_test = features_test[key]

            all_fields = set(f_ref.keys()) | set(f_test.keys())
            for field in all_fields:
                val_ref = f_ref.get(field, '<missing>')
                val_test = f_test.get(field, '<missing>')
                if val_ref != val_test:
                    differences[field].append({
                        'type': f_ref['type'],
                        'location': key,
                        'val_ref': val_ref,
                        'val_test': val_test,
                        gff_ref_path.name: val_ref,
                        gff_test_path.name: val_test
                    })

        if differences:
            for field, diffs in sorted(differences.items()):
                print(f"\nField: {field}")
                print(f"  Total differences: {len(diffs)}")
                for diff in diffs[:3]:
                    #if diff['type'] == 'CDS':
                    #    continue
                    print(f"    Type: {diff['type']}")
                    print(f"      Location: {diff['location']}")
                    print(f"      {gff_ref_path.name}: {diff[gff_ref_path.name]}")
                    print(f"      {gff_test_path.name}: {diff[gff_test_path.name]}")
                    print()
        
        else:
            print("\nNo differences found in common features!")


def build_parser():
    parser = argparse.ArgumentParser(
        description="Compare two GFF3 files and report differences (excluding ID/locus_tag)."
    )
    parser.add_argument("gff_ref", type=Path, help="Path to reference GFF3 file")
    parser.add_argument("gff_test", type=Path, help="Path to test GFF3 file")
    return parser


def main():
    args = build_parser().parse_args()

    if not args.gff_ref.exists():
        raise FileNotFoundError(f"Reference file not found: {args.gff_ref}")
    if not args.gff_test.exists():
        raise FileNotFoundError(f"Test file not found: {args.gff_test}")

    compare_gff3_files(args.gff_ref, args.gff_test)


if __name__ == "__main__":
    main()

