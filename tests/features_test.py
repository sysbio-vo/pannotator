#!/usr/bin/env python3

import argparse
from collections import defaultdict
from pathlib import Path
import csv
import sys

from BCBio import GFF


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
    """
    Print a full comparison of two GFF3 files (counts, ref-only/test-only
    features, and field-by-field diffs on common features).

    Returns True if the files are identical (no ref-only, test-only, or
    field differences), False otherwise -- used by batch mode to build the
    summary table without re-parsing or duplicating the diff logic.
    """
    features_ref = parse_gff3_to_dict(gff_ref_path)
    features_test = parse_gff3_to_dict(gff_test_path)

    keys_ref = features_ref.keys()
    keys_test = features_test.keys()

    only_in_ref = keys_ref - keys_test
    only_in_test = keys_test - keys_ref

    print(f"Total features in {gff_ref_path.name}: {len(keys_ref)}")
    print(f"Total features in {gff_test_path.name}: {len(keys_test)}")
    print(f"Common features: {len(keys_ref & keys_test)}")
    print('-'*80)

    # features only in ref file
    if only_in_ref:
        print(f"Only in {gff_ref_path.name}: {len(only_in_ref)}")
        print(f"Features:\n")
        for key in sorted(only_in_ref):
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
    if only_in_test:
        print(f"Only in {gff_test_path.name}: {len(only_in_test)}")
        print("Features:\n")
        for key in sorted(only_in_test):
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
    has_field_differences = False
    if common_keys:
        print("Difference in common features:\n")
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
            has_field_differences = True
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
            print("No differences found in common features!")

    return not (only_in_ref or only_in_test or has_field_differences)

def read_pairs(pairs_file: Path) -> list[tuple[Path, Path]]:
    """
    Read (ref_path, test_path) pairs from a TSV mapping file. Returns a
    list of (Path, Path) tuples. Raises ValueError on malformed rows.
    """
    with open(pairs_file, "r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        rows = list(reader)

    pairs = []
    for i, row in enumerate(rows, start=1):
        row = [cell for cell in row if cell != ""]  # tolerate trailing empty cells
        if not row:
            continue  # skip blank lines
        if len(row) != 2:
            raise ValueError(f"Row {i} in {pairs_file}: expected 2 columns, got {len(row)}: {row}")
        pairs.append((Path(row[0].strip()), Path(row[1].strip())))

    return pairs


def run_multiple_pairwise(input_file: Path) -> int:
    """
    Run compare_gff3_files() for every pair in a TSV mapping file, then
    print a summary table. Returns a process exit code (0 = all identical,
    1 = at least one difference or a missing path).
    """
    pairs = read_pairs(input_file)
    if not pairs:
        print(f"No pairs found in {input_file}.")
        return 1

    missing = []
    for ref_path, test_path in pairs:
        if not ref_path.exists():
            missing.append(str(ref_path))
        if not test_path.exists():
            missing.append(str(test_path))

    if missing:
        print("ERROR: the following paths do not exist:")
        for path in missing:
            print(f"  {path}")
        return 1

    results = []
    for ref_path, test_path in pairs:
        label = ref_path.name if ref_path.name == test_path.name else f"{ref_path.name} vs {test_path.name}"
        print(f"\n{'=' * 80}")
        print(label)
        print('=' * 80)
        identical = compare_gff3_files(ref_path, test_path)
        results.append((label, identical))

    print(f"\n{'=' * 80}")
    print("SUMMARY")
    print('=' * 80)
    print(f"{'Assembly':<40} {'Status':<8}")
    for label, identical in results:
        print(f"{label:<40} {'OK' if identical else 'DIFF':<8}")

    n_ok = sum(1 for _, identical in results if identical)
    print(f"\n{n_ok}/{len(results)} annotations identical.")

    return 0 if n_ok == len(results) else 1


def build_parser():
    parser = argparse.ArgumentParser(
        description="Compare GFF3 files and report differences (excluding ID/locus_tag)."
    )

    single = parser.add_argument_group(
        "single pair comparison",
        "Compare exactly one reference/test GFF3 pair supplied on the command line.",
    )
    single.add_argument("gff_ref", type=Path, nargs="?", default=None, help="Path to reference GFF3 file")
    single.add_argument("gff_test", type=Path, nargs="?", default=None, help="Path to test GFF3 file")

    multi = parser.add_argument_group(
        "multiple pairwise comparisons",
        "Compare multiple ref/test pairs listed in a TSV mapping file, instead of a single pair.",
    )
    multi.add_argument(
        "--input-file",
        type=Path,
        default=None,
        help="TSV file with no header and two columns per row: ref_path, test_path",
    )
    return parser


def main():
    parser = build_parser()
    args = build_parser().parse_args()

    # Validate input
    single_mode = args.gff_ref is not None or args.gff_test is not None
    multi_mode = args.input_file is not None

    if single_mode and multi_mode:
        parser.error("provide either <gff_ref> <gff_test> OR --input-file, not both")
    if not single_mode and not multi_mode:
        parser.error("provide either <gff_ref> <gff_test> OR --input-file")
    if single_mode and (args.gff_ref is None or args.gff_test is None):
        parser.error("both gff_ref and gff_test are required together")

    if multi_mode:
        exit_code = run_multiple_pairwise(args.input_file)
        sys.exit(exit_code)

    if not args.gff_ref.exists():
        raise FileNotFoundError(f"Reference file not found: {args.gff_ref}")
    if not args.gff_test.exists():
        raise FileNotFoundError(f"Test file not found: {args.gff_test}")

    compare_gff3_files(args.gff_ref, args.gff_test)


if __name__ == "__main__":
    main()
