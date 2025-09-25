import os
import pickle
import argparse

import utils as ut


def generate_pangenome(infile: str, outfile: str) -> None:
    clusters = ut.read_clusters(infile)

    # TODO: fill in
    pangenome = None

    ut.write_pickle(pangenome, outfile)


if __name__ == "__main__":
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description="Generate a pangenome index using homologous sequence clustering result.")
    parser.add_argument('-i', '--input', required=True, type=str, help='Clustering result as produced by MMseqs2')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output pangenome index')

    args = parser.parse_args()
