""" cluster_filter.py

This tool exports the clusters which contain the spectrun in specified project list cluster_id list.

Usage:
  cluster_filter.py --input=<results.clustering> --output=<filtered_results.clustering> 
                       [--project_ids=<ids.txt>]
                       [--min_size=<size>] [--max_size=<size>] [--min_ratio=<ratio>] [--max_ratio=<ratio>]
                       [--min_identified=<spectra>] [--max_identified=<spectra>]
  cluster_filter.py --cluster_ids=<ids.txt> --input=<results.clustering> --output=<filtered_results.clustering>
  cluster_filter.py (--help | --version)

Options:
  -i, --input=<clustering file>        Path to the .clustering result file to process.
  -o, --output=<features.txt>          Path to the output file that should be created. The output will
                                       be formatted as a tab-delimited text file.

  --project_ids=<prj_ids.txt>          If this parameter is set, the project ids are read from the
                                       specified file (one id per line) and only these clusters will
                                       be exported. All other filtering parameters are ignored if this
                                       parameter is set.
  --cluster_ids=<cluster_ids.txt>      If this parameter is set, the cluster ids are read from the
                                       specified file (one id per line) and only these clusters will
                                       be exported. All other filtering parameters are ignored if this
                                       parameter is set.
  --min_size=<size>                    The minimum size of a cluster to be reported.
  --max_size=<site>                    The maximum size of a cluster to be reported.
  --min_ratio=<ratio>                  The minimum ratio a cluster must have to be reported.
  --max_ratio=<ratio>                  The maximum ratio a cluster must have to be reported.
  --min_identified=<spectra>           May specify the minimum number of identified spectra a cluster must have.
  --max_identified=<spectra>           May specify the maximum number of identified spectra a cluster must have.
  -h, --help                           Print this help message.
  -v, --version                        Print the current version.
"""

import sys
import os
from docopt import docopt

# make the spectra_cluster packages available
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".."+ os.path.sep + "..")

from spectra_cluster.analyser.cluster_filter_analyser import ClusterFilter
import spectra_cluster.clustering_parser as clustering_parser



def create_analyser(arguments, output_file, set_params=True):
    """
    Creates an Cluster Filter analyser based on the command line
    parameters.
    :param arguments: The command line parameters
    :param output_file: File object opened to write to the output file
                        location
    :param set_params: If set the analysers parameters such as min / max
                       cluster size are set.
    :return: A Cluster Filter object
    """
    analyser = ClusterFilter(output_file)

    if not set_params:
        return analyser

    if arguments["--min_size"] is not None:
        analyser.min_size = int(arguments.get("--min_size"))
    if arguments["--max_size"] is not None:
        analyser.max_size = int(arguments.get("--max_size"))
    if arguments["--min_ratio"] is not None:
        analyser.min_ratio = float(arguments.get("--min_ratio"))
    if arguments["--max_ratio"] is not None:
        analyser.max_ratio = float(arguments.get("--max_ratio"))
    if arguments["--min_identified"] is not None:
        analyser.min_identified_spectra = int(arguments.get("--min_identified"))
    if arguments["--max_identified"] is not None:
        analyser.max_identified_spectra = int(arguments.get("--max_identified"))

    return analyser


def main():
    """
    Primary entry function for the CLI.
    :return:
    """
    arguments = docopt(__doc__, version='cluster_filter.py 1.0 BETA')

    # make sure the input file exists
    if not os.path.isfile(arguments['--input']):
        print("Error: Cannot find input file '" + arguments["--input"] + "'")
        sys.exit(1)

    # make sure the output file does not exist
    if os.path.isfile(arguments["--output"]):
        print("Error: Output file exists '" + arguments["--output"] + "'")
        sys.exit(1)

    project_ids = None
    cluster_ids = None

    if "--project_ids" in arguments and arguments["--project_ids"] is not None:
        project_ids = list()

        with open(arguments["--project_ids"], "r") as IN:
            for line in IN:
                project_ids.append(line.strip())
    if "--cluster_ids" in arguments and arguments["--cluster_ids"] is not None:
        cluster_ids = list()

        with open(arguments["--cluster_ids"], "r") as IN:
            for line in IN:
                cluster_ids.append(line.strip())


    with open(arguments["--output"] + ".part", "w") as OUT:
        # create the id transferer based on the settings
        analyser = create_analyser(arguments, OUT)

        # process all clusters
        parser = clustering_parser.ClusteringParser(arguments["--input"])

        print("Parsing input .clustering file...")
        for cluster in parser:
            # filter based on cluster ids if set
            if cluster_ids is not None:
                if cluster.id not in cluster_ids:
                    continue
            if project_ids is not None:
                analyser.remove_spectra_by_project(cluster,project_ids)
            analyser.process_cluster(cluster)
    os.rename(arguments["--output"] + ".part", arguments["--output"])
    print("Results written to " + arguments["--output"])


if __name__ == "__main__":
    main()
