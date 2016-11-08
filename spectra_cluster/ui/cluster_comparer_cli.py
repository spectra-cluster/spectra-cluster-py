"""cluster_comparer_cli

Command line interface to the spectra-cluster cluster as comparer tool. This tool compare two
cluster lists and give the statistics between them.

Usage:
  cluster_comparer_cli.py --input <results1.clustering,results2.clustering>
			  --output <comparer.txt>
                       [--min_size <size>] [--min_ratio <ratio>]
                       [--min_identified <spectra>]
                       [(--only_identified | --only_unidentified)]
  cluster_comparer_cli.py (--help | --version)

Options:
  -i, --input=<two clustering files>   Path to the two .clustering result files to process.
  -o, --output=<comparer.txt>          Path to the output file that should be created. The output will
                                       be formatted as a tab-delimited text file.
  --min_size=<size>                    The minimum size of a cluster to be reported. [default: 2]
  --min_ratio=<ratio>                  The minimum ratio a cluster must have to be reported.
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

import spectra_cluster.analyser.cluster_comparer as cluster_comparer 
import spectra_cluster.clustering_parser as clustering_parser


def create_analyser(arguments):
    """
    Creates an comparer analyser based on the command line
    parameters.
    :param arguments: The command line parameters
    :return: An Comparer object
    """
    analyser = cluster_comparer.ClusterListsComparer()

    if arguments["--only_identified"]:
        analyser.add_to_unidentified = False
    if arguments["--only_unidentified"]:
        analyser.add_to_identified = False

    if arguments["--min_size"]:
        analyser.min_size = int(arguments["--min_size"])
    if arguments["--min_ratio"]:
        analyser.min_ratio = float(arguments["--min_ratio"])

    if arguments["--min_identified"] is not None:
        analyser.min_identified_spectra = int(arguments["--min_identified"])


    return analyser


def write_results(identification_references, output_filename):
    """
    Writes the statistics as a tab delimited text file
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
    arguments = docopt(__doc__, version='cluster_comparer_cli 1.0 BETA')
#    print(arguments)
#    sys.exit(1)

    # make sure the input file exists
    input_files = arguments['--input'].split(',', 1)
    for input_file in input_files:
        if not os.path.isfile(input_file):      
            print("Error: Cannot find input file '" + input_file + "'")
            sys.exit(1)
    # make sure the output file does not exist
    if os.path.isfile(arguments["--output"]):
        print("Error: Output file exists '" + arguments["--output"] + "'")
        sys.exit(1)

    # create the id transferer based on the settings
    analyser = create_analyser(arguments)

    # process all clusters
    parser0 = clustering_parser.ClusteringParser(input_files[0])
    parser1 = clustering_parser.ClusteringParser(input_files[1])


    print("Parsing input .clustering file...")
    analyser.dealing_list = 0 #store the cluster in list0 
    for cluster in parser0:
        print("checking" + cluster.id)
        analyser.process_cluster(cluster)
    analyser.dealing_list = 1 #store the cluster in list1 
    for cluster in parser1:
        analyser.process_cluster(cluster)

    #do the compare
    analyser.compare() 

    #for debuging
    analyser.output_debug_info()

    # create the output file
#    write_results(analyser.identification_references, arguments["--output"])

    print("Results written to " + arguments["--output"])


if __name__ == "__main__":
    main()
