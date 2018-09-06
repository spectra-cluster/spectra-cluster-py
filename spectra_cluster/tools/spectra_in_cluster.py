"""
This tool creates a simple table (tab-delimited file) that contains the cluster ids in
the first column and the spectra ids in the second one.

Usage:
    spectra_in_cluster.py --input=<result.clustering> --output=<spectra_in_cluster.tsv>
    spectra_in_cluster.py (--help | --version)

Options:
    -i, --input=<result.clustering>        The clustering result file to process.
    -o, --output=<spectra_in_cluster.tsv>  Path to the result table.
    -h, --help                            Print this help message.
    -v, --version                         Print the current version.
"""

from docopt import docopt
import sys
import os

sys.path.insert(0, os.path.abspath('..') + os.path.sep + "..")

from spectra_cluster import clustering_parser


def main():
    arguments = docopt(__doc__, version='spectra_in_cluster.py 1.0 BETA')
    clustering_file = arguments["--input"]
    output_file = arguments["--output"]

    if not os.path.isfile(clustering_file):
        print("Error: Failed to find clustering result file '" + clustering_file + "'")
        sys.exit(1)

    if os.path.isfile(output_file):
        print("Error: Output file '" + output_file + "' exists.")
        sys.exit(1)

    parser = clustering_parser.ClusteringParser(clustering_file)

    print("Processing " + clustering_file + "...")

    with open(output_file, "w") as OUT:
        OUT.write("cluster_id\tspectrum_title\n")

        for cluster in parser:
            for spec_ref in cluster.get_spectra():
                OUT.write(cluster.id + "\t" + spec_ref.get_title() + "\n")

    print("Results written to '" + output_file + "'.")


if __name__ == "__main__":
    main()
