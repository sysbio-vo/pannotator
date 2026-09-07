# Description

The test suite shipped with Pannotator is designed to assess annotations of entire genomic datasets by comparing them to existing baseline (ground truth) data.

# How to use the test suite

First, you need to collect all reference (ground truth) GFF annotation files into a dedicated directory. Similarly, gather all annotations produced by Pannotator into a separate folder. Make sure that the annotations for the same samples have identical file names in both folders.

## Create an environment

To parse GFF files, the suite utilizes the `BCBio` package. Thus, you need to set up and activate a virtual environment before running the actual test.

```bash
python -m venv .pannotator_test_suite_venv
source .pannotator_test_suite_venv/bin/activate
pip install -r requirements.txt
```

## Run the test suite

To process the full dataset and print the feature similarity summary, simply run the test script.

```bash
python test_cds_search.py /path/to/reference/annotations /path/to/test/annotations > compare.log 2> compare.err
```

For a richer per-sample summary, append `verbose` as the third positional argument in the command. This will extend the logs by printing similarity statistics for each file pair individually.

```bash
python test_cds_search.py /path/to/reference/annotations /path/to/test/annotations verbose > compare.log 2> compare.err
```

In either case, the test suite generates two JSON files with the annotation matching summaries:

- `per_file_summary.json`
- `dataset_summary.json`

You can use those files to explore the GFF comparison results on your own.
