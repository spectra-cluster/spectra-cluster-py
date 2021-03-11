"""cluster_features_cli

Command line interface to the spectra-cluster cluster as features tool. This tool extracts the
clusters and samples found in the dataset and exports a table where each sample is represented
as a column and each cluster as a row. The cells then contain the number of spectra observed
from the given sample in the given cluster.

Usage:
  cluster_features_cli.py --input=<results.clustering> --output=<features.txt>
                       [--min_size=<size>] [--min_ratio=<ratio>]
                       [--min_identified=<spectra>] [--max_identified=<spectra>] 
                       [--output_matrix]
  cluster_features_cli.py (--help | --version)

Options:
  -i, --input=<clustering file>        Path to the .clustering result file to process.
  -o, --output=<features.txt>          Path to the output file that should be created. The output will
                                       be formatted as a tab-delimited text file.
  --min_size=<size>                    The minimum size of a cluster to be reported [default: 0]
  --min_ratio=<ratio>                  The minimum ratio a cluster must have to be reported [default: 0]
  --min_identified=<spectra>           May specify the minimum number of identified spectra a cluster must have [default: 0]
  --max_identified=<spectra>           May specify the maximum number of identified spectra a cluster may have [default: 1000000000]
  --output_matrix                      If set, the resulting file will only contain a matrix, with the first column
                                       as the cluster id, the second column the sample id, and the third column the
                                       number of spectra.
  -h, --help                           Print this help message.
  -v, --version                        Print the current version.
"""

import sys
import os
from docopt import docopt

# make the spectra_cluster packages available
sys.path.insert(0, os.path.abspath('..') + os.path.sep + "..")

from spectra_cluster.analyser.cluster_features import ClusterAsFeatures
import spectra_cluster.clustering_parser as clustering_parser


def create_analyser(arguments: dict, output_file):
    """
    Creates an ClusterAsFeatures analyser based on the command line
    parameters.
    :param arguments: The command line parameters
    :param output_file: File object opened to write to the output file
                        location
    :return: A ClusterAsFeatures object
    """
    analyser = ClusterAsFeatures(output_file)

    analyser.min_size = int(arguments.get("--min_size", 0))
    analyser.min_ratio = float(arguments.get("--min_ratio", 0))
    analyser.min_identified_spectra = int(arguments.get("--min_identified", 0))
    analyser.max_identified_spectra = int(arguments.get("--max_identified", 1000000000))

    return analyser


def main():
    """
    Primary entry function for the CLI.
    :return:
    """
    arguments = docopt(__doc__, version='id_transferer_cli 1.0 BETA')

    # make sure the input file exists
    if not os.path.isfile(arguments['--input']):
        print("Error: Cannot find input file '" + arguments["--input"] + "'")
        sys.exit(1)

    # make sure the output file does not exist
    if os.path.isfile(arguments["--output"]):
        print("Error: Output file exists '" + arguments["--output"] + "'")
        sys.exit(1)

    with open(arguments["--output"], "w") as OUT:
        # create the id transferer based on the settings
        analyser = create_analyser(arguments, OUT)

        # process all clusters
        parser = clustering_parser.ClusteringParser(arguments["--input"])

        print("Parsing input .clustering file...", end="", flush=True)
        processed_clusters = 0
        for cluster in parser:
            analyser.process_cluster(cluster)

            processed_clusters += 1
            if processed_clusters == 1000:
                print(".", end="", flush=True)
                processed_clusters = 0

    # add the header to the output file
    print("Adding header line...")
    analyser.add_resultfile_header(arguments["--output"])

    print("Results written to " + arguments["--output"])


if __name__ == "__main__":
    main()
