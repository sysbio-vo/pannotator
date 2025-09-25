import pickle

from typing import Any

# TODO: alternatively use sqlite for the pangenome index (potentially faster and more robust)
# NOTE: pickle loading might be slow

def read_pickle(file: str) -> dict:
    with open(file, 'rb') as f:
        res_dict = pickle.load(f)
    return res_dict

def write_pickle(obj: Any, outfile: str):
    with open(outfile, 'wb') as f:
        pickle.dump(obj, f)

def read_mmseqs_clusters(out_tsv: str) -> dict:
    with open(out_tsv, 'r') as f:
        clusters = f.read()
    return clusters

def read_clusters(file: str) -> dict:
    return read_mmseqs_clusters(file)
