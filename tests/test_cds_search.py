

import os
import sys
from pathlib import Path
import unittest
from BCBio import GFF
from concurrent.futures import ThreadPoolExecutor as Pool


def gff3_to_dict(gff_path: str, limit_info: dict | None = None) -> dict:
    res = dict()
    
    # limit_info = dict(gff_type=["CDS"]) # example

    with open(gff_path, 'r') as in_handle:
        for contig_rec in GFF.parse(in_handle, limit_info=limit_info):
            for contig_feature in contig_rec.features:

                res[contig_feature.id] = {
                    'type': contig_feature.type,
                    'location': str(contig_feature.location),
                    'qualifiers': contig_feature.qualifiers
                }

    return res


class TestCDSSearch(unittest.TestCase):
    REFERENCE_CDS_DIR = 'reference_data'
    TEST_CDS_DIR = 'test_data' # default values

    def process_test_file(self, filename: str):
        reference_path = os.path.join(self.REFERENCE_CDS_DIR, filename)
        test_path = os.path.join(self.TEST_CDS_DIR, filename)

        reference_dict = gff3_to_dict(reference_path, limit_info={'gff_type': ['CDS']})
        test_dict = gff3_to_dict(test_path, limit_info={'gff_type': ['CDS']})

        return reference_dict == test_dict


    def test_cds(self):
        # the files in the directories must have the same names
        test_files = sorted(os.listdir(self.TEST_CDS_DIR))

        with Pool() as pool:
            results = pool.map(self.process_test_file, test_files)

        self.assertTrue(all(results))       
        

if __name__ == '__main__':
    if len(sys.argv) > 1:
        print(sys.argv[1], sys.argv[2])
        TestCDSSearch.REFERENCE_CDS_DIR = sys.argv[1]
        TestCDSSearch.TEST_CDS_DIR = sys.argv[2]

        sys.argv = [sys.argv[0]]

    unittest.main()

