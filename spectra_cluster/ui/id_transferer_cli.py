"""id_transferer_cli

Command line interface to the spectra-cluster id_transferer tool. This tool can be
used to transfer identifications to other incorrectly or unidentified spectra.

Quality criteria can be specified in terms of minimum cluster size and minimum ratio (the
proportion of spectra identified as the most common peptide within a cluster). Additionally,
it is possible to specify the minimum number of identified spectra a cluster must have.

The tool then writes out all spectra references from the clusters that match the criteria
together with the most common peptide identification (PSM) of the respective cluster.

Usage:
  id_transferer_cli.py --input=<results.clustering> --output=<identifications.txt>
                       [--min_size=<size>] [--min_ratio=<ratio>]
                       [--min_identified=<spectra>]
                       [(--only_identified | --only_unidentified)]
  id_transferer_cli.py (--help | --version)

Options:
  -i, --input=<clustering file>        Path to the .clustering result file to process.
  -o, --output=<identifications.txt>   Path to the output file that should be created. The output will
                                       be formatted as a tab-delimited text file.
  --min_size=<size>                    The minimum size of a cluster to be reported. [default: 5]
  --min_ratio=<ratio>                  The minimum ratio a cluster must have to be reported. [default: 0.7]
  --min_identified=<spectra>           May specify the minimum number of identified spectra a cluster must have.
  --only_identified                    If set, only identified spectra will be reported.
  --only_unidentified                  If set, only unidentified spectra will be reported.
  -h, --help                           Print this help message.
  -v, --version                        Print the current version.
"""

import sys
import os
from docopt import docopt

# make the spectra_cluster packages available
sys.path.insert(0, os.path.abspath('..') + os.path.sep + "..")

import spectra_cluster.analyser.id_transferer as id_transferer
import spectra_cluster.clustering_parser as clustering_parser


def create_analyser(arguments):
    """
    Creates an IdTransferer analyser based on the command line
    parameters.
    :param arguments: The command line parameters
    :return: An IdTransferer object
    """
    analyser = id_transferer.IdTransferer()

    if arguments["--only_identified"]:
        analyser.add_to_unidentified = False
    if arguments["--only_unidentified"]:
        analyser.add_to_identified = False

    analyser.min_size = int(arguments["--min_size"])
    analyser.min_ratio = float(arguments["--min_ratio"])

    if arguments["--min_identified"] is not None:
        analyser.min_identified_spectra = int(arguments["--min_identified"])

    return analyser


def write_results(identification_references, output_filename):
    """
    Writes the identification references as a tab delimited text file
    to the specified path.
    :param identification_references: List of identification references.
    :param output_filename: Path to the output filename
    :return:
    """
    with open(output_filename, "w") as writer:
        writer.write("filename\tspec_id\tsequence\n")

        for id_ref in identification_references:
            psm_string = ";".join([str(p) for p in identification_references.psms])

            writer.write(id_ref.filename + "\t" + id_ref.spec_id + "\t" + psm_string + "\n")


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

    # create the id transferer based on the settings
    analyser = create_analyser(arguments)

    # process all clusters
    parser = clustering_parser.ClusteringParser(arguments["--input"])

    print("Parsing input .clustering file...")
    for cluster in parser:
        analyser.process_cluster(cluster)

    # create the output file
    write_results(analyser.identification_references, arguments["--output"])

    print("Results written to " + arguments["--output"])


if __name__ == "__main__":
    main()
