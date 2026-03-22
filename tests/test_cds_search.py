import os
import re
import sys
import unittest
from concurrent.futures import ThreadPoolExecutor as Pool

from BCBio import GFF

CONTIG_ID = r"(?P<contig_id>\d+)$"


def gff3_to_dict(gff_path: str, limit_info: dict | None = None) -> dict:
    res = dict()

    # limit_info = dict(gff_type=["CDS"]) # example

    with open(gff_path, "r") as in_handle:
        for contig_rec in GFF.parse(in_handle, limit_info=limit_info):
            contig_id = contig_rec.id
            for contig_feature in contig_rec.features:

                res[contig_feature.id] = {
                    "contig_id": contig_id,
                    "type": contig_feature.type,
                    "location": str(contig_feature.location),
                    "qualifiers": contig_feature.qualifiers,
                }

    return res


# elaborate helpers to estimate prediction accuracy


def use_contig_id_and_location_as_locus_tag(gff_dict: dict) -> dict:
    res = dict()
    for k, v in gff_dict.items():
        contig_num = re.search(CONTIG_ID, v["contig_id"])
        if contig_num is not None:
            contig_num = int(contig_num.group("contig_id"))
        res[f"{v['location']}_{contig_num}"] = v
    return res


def calculate_dict_keys_similarity(ref_dict, test_dict, verbose=False):
    ref_keys = set(ref_dict.keys())
    test_keys = set(test_dict.keys())

    recall_pct = len(ref_keys.intersection(test_keys)) / len(ref_keys) * 100

    if verbose:
        print(
            f"""
Total number of reference features: {len(ref_keys)}
Total number of features found with Pannotator: {len(test_keys)}
Overlap: {len(ref_keys.intersection(test_keys))} ({recall_pct:.2f}% of reference features)
Number of found features missing in reference set: {len(test_keys - ref_keys)}
Number of reference features that weren't found: {len(ref_keys - test_keys)}
"""
        )

    return recall_pct


def convert_incorrectly_parsed_list_to_str(features_list):
    if len(features_list) > 1:
        feature_str = ",".join(features_list)
    else:
        feature_str = features_list[0]
    return feature_str


def compare_incorrectly_parsed_list_features(ref_feature, test_feature):
    assert isinstance(ref_feature, list) and len(ref_feature) > 0
    assert isinstance(test_feature, list) and len(ref_feature) > 0

    ref_feature_str = convert_incorrectly_parsed_list_to_str(ref_feature)
    test_feature_str = convert_incorrectly_parsed_list_to_str(test_feature)

    return test_feature_str == ref_feature_str


def calcule_important_field_similarity(ref_dict, test_dict):

    type_errors = []
    name_errors = []
    dbxref_errors = []

    for k in test_dict.keys():
        if k not in ref_dict.keys():
            continue
        test_feature = test_dict[k]
        ref_feature = ref_dict[k]

        # important fields: type, qualifiers.Name, qualifiers.product, qualifieres.Dbxref
        type_eq = test_feature["type"] == ref_feature["type"]
        if not type_eq:
            type_errors.append((k, test_feature["type"], ref_feature["type"]))

        if "Name" in test_feature["qualifiers"] and "Name" in ref_feature["qualifiers"]:
            name_eq = compare_incorrectly_parsed_list_features(
                test_feature["qualifiers"]["Name"], ref_feature["qualifiers"]["Name"]
            )
            if not name_eq:
                type_errors.append(
                    (
                        k,
                        convert_incorrectly_parsed_list_to_str(test_feature["qualifiers"]["Name"]),
                        convert_incorrectly_parsed_list_to_str(ref_feature["qualifiers"]["Name"]),
                    )
                )

        if "Dbxref" in test_feature["qualifiers"] and "Dbxref" in ref_feature["qualifiers"]:
            dbxref_eq = sorted(test_feature["qualifiers"]["Dbxref"]) == sorted(ref_feature["qualifiers"]["Dbxref"])
            if not dbxref_eq:
                dbxref_errors.append(
                    (k, sorted(test_feature["qualifiers"]["Dbxref"]), sorted(ref_feature["qualifiers"]["Dbxref"]))
                )

    return type_errors, name_errors, dbxref_errors


def calculate_overall_similarity(ref_dict, test_dict, sample_name: str = ""):
    mismatches = dict()
    total_num_of_field_matches = 0
    total_num_of_fields = 0
    for k, v in test_dict.items():
        if not isinstance(v, dict):
            if k in ref_dict.keys():
                if test_dict[k] != ref_dict[k]:
                    if f"{sample_name}_{k}" in mismatches.keys():
                        mismatches[f"{sample_name}_{k}"].append((test_dict[k], ref_dict[k]))
                    else:
                        mismatches[f"{sample_name}_{k}"] = [(test_dict[k], ref_dict[k])]
                else:
                    total_num_of_field_matches += test_dict[k] == ref_dict[k]
            total_num_of_fields += 1
        elif k in ref_dict.keys() and isinstance(ref_dict[k], dict):
            subdict_matches, subdict_fields, subdict_mismatches = calculate_overall_similarity(
                ref_dict[k], test_dict[k], sample_name
            )
            total_num_of_field_matches += subdict_matches
            total_num_of_fields += subdict_fields
            for miss_k, miss_v in subdict_mismatches.items():
                if miss_k in mismatches.keys():
                    mismatches[miss_k] += miss_v
                else:
                    mismatches[miss_k] = miss_v
        else:
            total_num_of_fields += 1

    return total_num_of_field_matches, total_num_of_fields, mismatches


def print_top_n_mismatches_of_each_type(mismatches_dict: dict, n: int = 5) -> None:
    for k, v in mismatches_dict.items():
        mismatches_print_str = "\n".join(
            map(lambda mismatch: f'Predicted value is "{mismatch[0]}", reference value is "{mismatch[1]}"', v[:n])
        )
        print(f"{n} first mismatches for field {k}: {mismatches_print_str} \n")


def get_print_str_for_top_n_important_field_mismatches(field_errors: list, field_name: str, n: int = 5) -> None:
    return "\n".join(
        list(
            map(
                lambda mismatch: f'Mismatch at location {mismatch[0]}: predicted {field_name} \
is "{mismatch[1]}", reference value is "{mismatch[2]}"',
                field_errors[:n],
            )
        )
    )


def calculate_cds_content_similarity(
    ref_dict, test_dict, mismatch_n_to_print: int = 5, verbose=False, sample_name: str = ""
):
    # calculate retained region sets similarity
    recall_pct = calculate_dict_keys_similarity(ref_dict, test_dict, verbose)

    # calculate important fields similarity
    type_errors, name_errors, dbxref_errors = calcule_important_field_similarity(ref_dict, test_dict)

    type_error_rate = len(type_errors) / len(test_dict)
    name_error_rate = len(name_errors) / len(test_dict)
    dbxref_error_rate = len(dbxref_errors) / len(test_dict)

    total_features_num = len(test_dict)

    if verbose:

        print(
            f"""
{f'Sample {sample_name}' if sample_name else ''}
Total number of sequence type mismatches: {len(type_errors)} ({type_error_rate * 100:.2f}%)
{get_print_str_for_top_n_important_field_mismatches(type_errors, 'type', mismatch_n_to_print)}

Total number of product name mismatches: {len(name_errors)} ({name_error_rate * 100:.2f}%)
{get_print_str_for_top_n_important_field_mismatches(name_errors, 'name', mismatch_n_to_print)}

Total number of Dbxref name mismatches: {len(dbxref_errors)} ({dbxref_error_rate * 100:.2f}%)
{get_print_str_for_top_n_important_field_mismatches(dbxref_errors, 'Dbxref', mismatch_n_to_print)}
"""
        )

    total_field_matches, total_field_num, mismatches = calculate_overall_similarity(ref_dict, test_dict, sample_name)

    if verbose:
        print(
            f"""
Overall result similarity: {total_field_matches / total_field_num * 100: .2f}%
"""
        )
        print_top_n_mismatches_of_each_type(mismatches, mismatch_n_to_print)

    important_field_mean_error_rate = sum([type_error_rate, name_error_rate, dbxref_error_rate]) / 3

    summary_dict = dict(
        recall_pct=recall_pct,
        type_errors=type_errors,
        name_errors=name_errors,
        dbxref_errors=dbxref_errors,
        important_field_mean_error_rate=important_field_mean_error_rate,
        type_error_rate=type_error_rate,
        name_error_rate=name_error_rate,
        dbxref_error_rate=dbxref_error_rate,
        total_field_matches=total_field_matches,
        total_field_num=total_field_num,
        mismatches=mismatches,
        total_features_num=total_features_num,
    )

    return summary_dict


def print_aggregated_summary(summary_dicts_list: list[dict], mismatch_n_to_print: int = 5) -> None:
    recall_pcts = []
    important_field_mean_errors = []

    total_field_matches_in_dataset = 0
    total_field_num_in_dataset = 0
    total_features_num_in_dataset = 0

    type_errors = []
    name_errors = []
    dbxref_errors = []

    type_error_rates = []
    name_error_rates = []
    dbxref_error_rates = []

    mismatches = dict()

    for summary_dict in summary_dicts_list:
        if summary_dict is None:
            continue
        recall_pcts.append(summary_dict["recall_pct"])
        important_field_mean_errors.append(summary_dict["important_field_mean_error_rate"])
        total_field_matches_in_dataset += summary_dict["total_field_matches"]
        total_field_num_in_dataset += summary_dict["total_field_num"]
        total_features_num_in_dataset += summary_dict["total_features_num"]

        type_errors += summary_dict["type_errors"]
        type_error_rates.append(len(summary_dict["type_errors"]) / summary_dict["total_features_num"])

        name_errors += summary_dict["name_errors"]
        name_error_rates.append(len(summary_dict["name_errors"]) / summary_dict["total_features_num"])

        dbxref_errors += summary_dict["dbxref_errors"]
        dbxref_error_rates.append(len(summary_dict["dbxref_errors"]) / summary_dict["total_features_num"])

        for k, v in summary_dict["mismatches"].items():
            if k in mismatches:
                mismatches[k] += v
            else:
                mismatches[k] = v

    type_error_rate = len(type_errors) / total_features_num_in_dataset
    name_error_rate = len(name_errors) / total_features_num_in_dataset
    dbxref_error_rate = len(dbxref_errors) / total_features_num_in_dataset

    max_important_field_mismatch_rate = max(max(type_error_rates), max(name_error_rates), max(dbxref_error_rates))

    print(
        f"""
Mean found features recall percentage: {sum(recall_pcts) / len(recall_pcts):.2f}%
Minimum features recall percentage: {min(recall_pcts):.2f}%
Total number of sequence type mismatches: {len(type_errors)} ({type_error_rate * 100:.2f}%)
{get_print_str_for_top_n_important_field_mismatches(type_errors, 'type', mismatch_n_to_print)}

Total number of product name mismatches: {len(name_errors)} ({name_error_rate * 100:.2f}%)
{get_print_str_for_top_n_important_field_mismatches(name_errors, 'name', mismatch_n_to_print)}

Total number of Dbxref mismatches: {len(dbxref_errors)} ({dbxref_error_rate * 100:.2f}%)
{get_print_str_for_top_n_important_field_mismatches(dbxref_errors, 'Dbxref', mismatch_n_to_print)}

Maximum important field mismatch rate: {max_important_field_mismatch_rate * 100:.2f}%

Overall result similarity: {total_field_matches_in_dataset / total_field_num_in_dataset * 100:.2f}%
"""
    )

    print_top_n_mismatches_of_each_type(mismatches, mismatch_n_to_print)

    return max_important_field_mismatch_rate


class TestCDSSearch(unittest.TestCase):
    REFERENCE_CDS_DIR = "reference_data"
    TEST_CDS_DIR = "test_data"  # default values
    VERBOSE = False
    ADMISSIBLE_IMPORTANT_FIELD_ERROR_RATE = 0.05

    def process_test_file(self, filename: str, mismatch_n_to_print: int = 5):
        reference_path = os.path.join(self.REFERENCE_CDS_DIR, filename)
        test_path = os.path.join(self.TEST_CDS_DIR, filename)

        if not os.path.exists(reference_path) or not os.path.exists(test_path):
            return None

        reference_dict = gff3_to_dict(reference_path, limit_info=None)
        test_dict = gff3_to_dict(test_path, limit_info=None)

        ref_dict_location_and_contig = use_contig_id_and_location_as_locus_tag(reference_dict)
        test_dict_location_and_contig = use_contig_id_and_location_as_locus_tag(test_dict)

        # WARNING: this will print summary for each file in the dataset
        # TODO: limit or aggregate the output
        summary_dict = calculate_cds_content_similarity(
            ref_dict_location_and_contig, test_dict_location_and_contig, mismatch_n_to_print, self.VERBOSE, filename
        )

        return summary_dict

    def test_cds(self):
        # the files in the directories must have the same names
        test_files = sorted(os.listdir(self.TEST_CDS_DIR))

        with Pool() as pool:
            sumary_dicts = pool.map(self.process_test_file, test_files)

        max_important_field_error_rate = print_aggregated_summary(sumary_dicts)

        self.assertTrue(max_important_field_error_rate <= self.ADMISSIBLE_IMPORTANT_FIELD_ERROR_RATE)


if __name__ == "__main__":
    if len(sys.argv) == 3:
        TestCDSSearch.REFERENCE_CDS_DIR = sys.argv[1]
        TestCDSSearch.TEST_CDS_DIR = sys.argv[2]

        sys.argv = [sys.argv[0]]

    elif len(sys.argv) >= 4:
        TestCDSSearch.REFERENCE_CDS_DIR = sys.argv[1]
        TestCDSSearch.TEST_CDS_DIR = sys.argv[2]

        verbose = sys.argv[3]
        if verbose.lower() == "verbose":
            TestCDSSearch.VERBOSE = True

        sys.argv = [sys.argv[0]]

    unittest.main()
