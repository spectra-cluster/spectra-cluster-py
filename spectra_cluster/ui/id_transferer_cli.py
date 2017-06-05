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
                       [--moff_compatible]
                       [--fasta=<species.fasta>]
                       [--min_size=<size>] [--min_ratio=<ratio>]
                       [--min_identified=<spectra>]
                       [(--only_identified | --only_unidentified)]
                       [--return_all_identified]
  id_transferer_cli.py (--help | --version)

Options:
  -i, --input=<clustering file>        Path to the .clustering result file to process.
  -o, --output=<identifications.txt>   Path to the output file that should be created. The output will
                                       be formatted as a tab-delimited text file.
  --moff_compatible                    If this parameter is specified, the output will be formatted
                                       according to the requirements of the "Modest Feature Finder" tool.
                                       This also requires protein inference to be performed which is the
                                       reason why a FASTA database must be supplied. Here, the smallest set
                                       of proteins required to explain the peptides will be reported.
  --fasta=<species.fasta>              If a fasta file is supplied, all peptides will be mapped to the respective
                                       proteins and an additional "protein" column will be added to the output
                                       file. MoFF compatible files will always contain the required protein column.
  --min_size=<size>                    The minimum size of a cluster to be reported. [default: 5]
  --min_ratio=<ratio>                  The minimum ratio a cluster must have to be reported. [default: 0.7]
  --min_identified=<spectra>           May specify the minimum number of identified spectra a cluster must have.
  --only_identified                    If set, only identified spectra will be reported.
  --only_unidentified                  If set, only unidentified spectra will be reported.
  --return_all_identified              If set, originally identified spectra are returned irrespective of in which
                                       cluster they are present. Additionally, if this option is set together with
                                       the --only_unidentified option, all originally identified spectra are returned
                                       unchanged.
  -h, --help                           Print this help message.
  -v, --version                        Print the current version.
"""

import sys
import os
from docopt import docopt
from maspy_resources import peptide_mapping

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
    analyser = id_transferer.IdTransferer(True, True, False)

    if arguments["--only_identified"]:
        analyser.add_to_unidentified = False
    if arguments["--only_unidentified"]:
        analyser.add_to_identified = False
    if arguments["--return_all_identified"]:
        analyser.include_all_identified = True

    analyser.min_size = int(arguments["--min_size"])
    analyser.min_ratio = float(arguments["--min_ratio"])

    if arguments["--min_identified"] is not None:
        analyser.min_identified_spectra = int(arguments["--min_identified"])

    return analyser


def write_results(identification_references, peptide_mappings, output_filename):
    """
    Writes the identification references as a tab delimited text file
    to the specified path.

    :param identification_references: List of identification references.
    :param peptide_mappings: Mappings of peptides to proteins (dict with peptide sequence as key)
    :param output_filename: Path to the output filename
    """
    with open(output_filename, "w") as writer:
        if peptide_mappings is not None:
            writer.write("filename\tspec_id\tsequence\tchanged_by_clustering\tprotein\n")
        else:
            writer.write("filename\tspec_id\tsequence\tchanged_by_clustering\n")

        for id_ref in identification_references:
            psm_string = ";".join([str(p) for p in id_ref.psms])

            changed_by_clustering = "false"
            if id_ref.changed_through_clustering:
                changed_by_clustering = "true"

            fields = [id_ref.filename, id_ref.spec_id, psm_string, changed_by_clustering]

            if peptide_mappings is not None:
                sequence = id_ref.psms[0].sequence
                proteins = peptide_mappings[sequence] if sequence in peptide_mappings else ""
                sorted_proteins = list(proteins)
                sorted_proteins.sort()
                fields.append(";".join(sorted_proteins))

            writer.write("\t".join(fields) + "\n")


def write_moff_results(identification_references, peptide_mappings, output_filename):
    """
    Writes the identification references as a tab delimited text file
    with compatible formatting required by the "Modest Feature Finder" tool.

    :param identification_references: List of identification references.
    :param peptide_mappings: Mappings of peptides to proteins (dict with peptide sequence as key)
    :param output_filename: Path to the output filename
    """
    with open(output_filename, "w") as writer:
        writer.write("peptide\tprot\tmod_peptide\trt\tmz\tmass\tcharge\n")

        for id_ref in identification_references:
            fields = list()

            # only report the first PSM
            psm = id_ref.psms[0]
            fields.append(psm.sequence)

            # protein mapping
            if peptide_mappings is not None and psm.sequence in peptide_mappings:
                proteins = list(peptide_mappings[psm.sequence])
                proteins.sort()

                fields.append(";".join(proteins))
            else:
                fields.append("")

            mod_peptide_string = str(psm)
            fields.append(mod_peptide_string)

            rt = id_ref.spectrum.get_property("RT")
            if rt is None:
                raise Exception("Missing retention time information for spectrum object. "
                                "Please cluster your data with the spectra-cluster-cli tool "
                                "version >= 1.0.4")

            fields.append(rt)

            fields.append(str(id_ref.spectrum.precursor_mz))

            fields.append(str(id_ref.spectrum.get_mass()))

            fields.append(str(id_ref.spectrum.charge))

            writer.write("\t".join(fields) + "\n")


def main():
    """
    Primary entry function for the CLI.

    """
    arguments = docopt(__doc__, version='id_transferer_cli 1.0 BETA')

    # make sure the input file exists
    input_file = arguments['--input']
    if not os.path.isfile(input_file):
        print("Error: Cannot find input file '" + input_file + "'")
        sys.exit(1)

    # make sure the output file does not exist
    output_file = arguments["--output"]
    if os.path.isfile(output_file):
        print("Error: Output file exists '" + output_file + "'")
        sys.exit(1)

    fasta_file = arguments["--fasta"]

    if fasta_file is not None:
        if not os.path.isfile(fasta_file):
            print("Error: Cannot find FASTA file '" + fasta_file + "'")
            sys.exit(1)

    # create the id transferer based on the settings
    analyser = create_analyser(arguments)

    # process all clusters
    parser = clustering_parser.ClusteringParser(input_file)

    print("Parsing input .clustering file...")
    for cluster in parser:
        analyser.process_cluster(cluster)

    # perform protein inference
    if fasta_file is not None:
        print("Doing protein inference...")
        all_peptides = set()
        for id_ref in analyser.identification_references:
            for psm in id_ref.psms:
                all_peptides.add(psm.sequence)

        peptide_mappings = peptide_mapping.the_magic_mapping_function(all_peptides, fasta_file)
    else:
        peptide_mappings = None

    # create the output file
    if arguments["--moff_compatible"]:
        write_moff_results(identification_references=analyser.identification_references,
                           peptide_mappings=peptide_mappings,
                           output_filename=output_file)
    else:
        write_results(identification_references=analyser.identification_references,
                      peptide_mappings=peptide_mappings,
                      output_filename=output_file)

    print("Results written to " + output_file)


if __name__ == "__main__":
    main()
