"""clustering_stats

Extracts basic statistics (ie. number of clusters, incorrectly clustered spectra) from .clustering files.
This script only creates meaningful results if the .clustering file contains identification data which
will be used to evaluate correctly and incorrectly clustered spectra.

Usage:
    clustering_stats --output <stats.tsv> --min_size <3> CLUSTERING_FILE...
    clustering_stats (--help | --usage)

Options:
  -o, --output <stats.tsv>          Name of the result file to be created.
  -s, --min_size <3>                The minimum size a cluster must have to be evaluated [default: 3]
  -h, --help                        Show this help message.
  --usage                           Show usage information
"""

import os
import sys
from docopt import docopt

# make the spectra_cluster packages available
sys.path.insert(0, os.path.abspath('..') + os.path.sep + "..")

import spectra_cluster.clustering_parser


class ClusteringStatistics:
    """
    This class is only intended to hold basic information about a .clustering file
    """
    def __init__(self, filename, min_size, total_clusters, total_spectra, identified_spectra,
                 clustered_spectra, correct_spectra, incorrect_spectra, min_size_clusters, clustered_identified_spectra):
        self.filename = filename
        self.min_size = min_size
        self.total_spectra = total_spectra
        self.total_clusters = total_clusters
        self.identified_spectra = identified_spectra
        self.clustered_spectra = clustered_spectra
        self.correct_spectra = correct_spectra
        self.incorrect_spectra = incorrect_spectra
        self.min_size_clusters = min_size_clusters
        self.clustered_identified_spectra = clustered_identified_spectra


def determine_clustering_accuracy(clustering_file, min_size=3):
    """
    Determines the clustering accuracy and number of clustered spectra
    in the analysed .clustering file

    :param clustering_file: The .clustering file to process
    :param min_size: The minimum sizes a cluster must have to be evaluated
    :return: dict(size => (rel. clustered spectra, accuracy))
    """
    parser = spectra_cluster.clustering_parser.ClusteringParser(clustering_file)

    total_spectra = 0
    clustered_spectra = 0
    clustered_identified_spectra = 0
    correct_spectra = 0
    incorrect_spectra = 0
    total_clusters = 0
    identified_spectra = 0
    min_size_clusters = 0

    # count correct and incorrect spectra per size
    for cluster in parser:
        total_spectra += cluster.n_spectra
        identified_spectra += cluster.identified_spectra
        total_clusters += 1

        if cluster.identified_spectra >= min_size:
            min_size_clusters += 1
            clustered_spectra += cluster.n_spectra
            clustered_identified_spectra += cluster.identified_spectra
            correct_spectra += cluster.max_il_ratio * cluster.identified_spectra
            incorrect_spectra += (1 - cluster.max_il_ratio) * cluster.identified_spectra

    return ClusteringStatistics(filename=clustering_file, min_size=min_size, total_clusters=total_clusters,
                                total_spectra=total_spectra, identified_spectra=identified_spectra,
                                clustered_spectra=clustered_spectra,
                                correct_spectra=correct_spectra, incorrect_spectra=incorrect_spectra,
                                min_size_clusters=min_size_clusters,
                                clustered_identified_spectra=clustered_identified_spectra)

def main():
    args = docopt(__doc__)

    result_file = args["--output"]
    min_size = int(args["--min_size"])

    if os.path.isfile(result_file):
        print("Error: Output file exists")
        sys.exit(1)

    # create the clustering statistics
    stats = list()
    for clustering_file in args["CLUSTERING_FILE"]:
        print("Processing " + clustering_file + "...")
        stats.append(determine_clustering_accuracy(clustering_file, min_size=min_size))

    # create the result file
    with open(result_file, "w") as writer:
        writer.write("filename\tmin_size\ttotal_clusters\ttotal_spectra\tidentified_spectra\t" +
                     "clustered_spectra\tclustered_identified_spectra\tcorrect_spectra\tincorrect_spectra\t" +
                     "min_sized_clusters\n")

        for result in stats:
            writer.write("\t".join([
                result.filename,
                str(result.min_size),
                str(result.total_clusters),
                str(result.total_spectra),
                str(result.identified_spectra),
                str(result.clustered_spectra),
                str(result.clustered_identified_spectra),
                str(result.correct_spectra),
                str(result.incorrect_spectra),
                str(result.min_size_clusters)
            ]) + "\n")

    print("Results written to " + result_file)


if __name__ == "__main__":
    main()
