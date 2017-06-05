"""cluster_parameter_extractor

Extracts basic parameters about the clusters found in a .clustering result file and
writes these parameters to a tab-delimited file.

Usage:
  cluster_parameter_extractor.py --input=<results.clustering> --output=<parameters.txt>
                                 [--synthetic_peptides]
  cluster_parameter_extractor.py (--help | --version)

Options:
  -i, --input=<clustering file>        Path to the .clustering result file to process.
  -o, --output=<parameters.txt>        Path to the output file that should be created.
  --synthetic_peptides                 If this option is specified, all spectra from the
                                       dataset on synthetic peptides (PXD004732) will be
                                       analysed separately.
  -h, --help                           Print this help message.
  -v, --version                        Print the current version.
"""

import sys
import os
from docopt import docopt

# make the spectra_cluster packages available
sys.path.insert(0, os.path.abspath('..') + os.path.sep + "..")

import spectra_cluster.clustering_parser as clustering_parser


def create_sequence_string(cluster):
    """
    Creates a string containing the cluster's sequences with their
    respective spectrum counts.

    :param cluster: The cluster to process
    :return: A string containing the cluster's sequences with their counts.
    """
    sequence_strings = list()

    for sequence in cluster.sequence_counts.keys():
        sequence_strings.append(sequence + ":" + str(cluster.sequence_counts.get(sequence)))

    final_string = "[" + ",".join(sequence_strings) + "]"

    return final_string


def process_cluster(cluster):
    """
    Extracts the basic properties from the cluster object and
    returns a tab-delimited string that matches the output format.

    :param cluster: Cluster object to extract the parameters from
    :return: A tab-delimited string matching the output format
    """
    result_fields = [cluster.id, str(cluster.precursor_mz), str(cluster.n_spectra),
                     str(cluster.identified_spectra), str(cluster.unidentified_spectra),
                     str(cluster.max_ratio), str(cluster.max_il_ratio)]

    # precursor m/z range
    spectra = cluster.get_spectra()
    min_mz = sys.float_info.__getattribute__("max")
    max_mz = 0

    for spectrum in spectra:
        if min_mz > spectrum.precursor_mz:
            min_mz = spectrum.precursor_mz
        if max_mz < spectrum.precursor_mz:
            max_mz = spectrum.precursor_mz

    mz_range = max_mz - min_mz

    result_fields.append(str(mz_range))

    # sequences
    sequence_string = create_sequence_string(cluster)
    result_fields.append(sequence_string)

    # max sequence
    if len(cluster.max_sequences) == 0:
        result_fields.append("NA")
        result_fields.append("NA")
        result_fields.append("NA")
    if len(cluster.max_sequences) > 0:
        result_fields.append(cluster.max_sequences[0])
        result_fields.append(str(cluster.sequence_counts.get(cluster.max_sequences[0])))
        # TODO: add modifications to output
        result_fields.append("NA")

    # second max sequence
    if len(cluster.sequence_counts) == 1:
        # cluster only contains 1 sequence
        result_fields.append("NA")
        result_fields.append("NA")
        result_fields.append("NA")
    elif len(cluster.max_sequences) > 1:
        # there is more than one max sequence
        result_fields.append(cluster.max_sequences[1])
        result_fields.append(str(cluster.sequence_counts.get(cluster.max_sequences[1])))
        # TODO: add modifications to output
        result_fields.append("NA")
    else:
        # identify the second most common sequence
        counts = list(cluster.sequence_counts.values())
        counts.remove(cluster.sequence_counts.get(cluster.max_sequences[0]))
        second_max_count = max(counts)

        sequence_found = False

        for sequence in cluster.sequence_counts.keys():
            if cluster.sequence_counts.get(sequence) == second_max_count:
                result_fields.append(sequence)
                sequence_found = True
                break

        # this should never happen
        if not sequence_found:
            raise Exception("Failed to identify second most common sequence for cluster " + cluster.id)

        result_fields.append(str(second_max_count))
        # TODO: add modifications
        result_fields.append("NA")

    # n_input_files
    input_files = set()

    for spectrum in spectra:
        spectrum_filename = spectrum.get_filename()

        if spectrum_filename is not None:
            input_files.add(spectrum_filename)

    if len(input_files) < 1:
        result_fields.append("NA")
    else:
        result_fields.append(str(len(input_files)))

    return "\t".join(result_fields)


def process_synthetic_peptides(cluster):
    """
    Extracts information on identification from the dataset
    PXD004732 which contains reliable identifications of
    synthetic peptides.

    :param cluster: The cluster to process
    :return: A tab-delimited string containing the number of
             synthetic peptide spectra in the cluster, their
             max ratio and the max sequence.
    """
    total_spectra = 0
    sequence_counts = dict()

    for spectrum in cluster.get_spectra():
        # ignore all none synthetic spectra
        if "PXD004732" not in spectrum.get_title():
            continue

        total_spectra += 1

        # spectra from this dataset only contain one psm
        if len(spectrum.psms) != 1:
            raise Exception("Unidentified spectrum from project PXD004732 encountered: " + spectrum.get_title())

        psm = [p for p in spectrum.psms][0]

        if psm.sequence not in sequence_counts:
            sequence_counts[psm.sequence] = 1
        else:
            sequence_counts[psm.sequence] += 1

    # done if no synthetic peptides were found
    if total_spectra == 0:
        return "0\tNA\tNA"

    # get the max count
    max_count = max(sequence_counts.values())

    # get the matching sequences
    max_sequences = [s for s in sequence_counts if sequence_counts[s] == max_count]

    # calculate the ratio
    max_ratio = max_count / total_spectra

    return str(total_spectra) + "\t" + str(max_ratio) + "\t" + ";".join(max_sequences)


def main():
    """
    Primary entry function for the CLI.
    :return:
    """
    arguments = docopt(__doc__, version='cluster_parameter_extractor 1.0 BETA')

    input_file = arguments['--input']
    output_file = arguments["--output"]
    process_synthetic = arguments["--synthetic_peptides"]

    # make sure the input file exists
    if not os.path.isfile(input_file):
        print("Error: Cannot find input file '" + input_file + "'")
        sys.exit(1)

    # make sure the output file does not exist
    if os.path.isfile(output_file):
        print("Error: Output file exists '" + output_file + "'")
        sys.exit(1)

    with open(output_file, "w") as OUT:
        # write the header
        OUT.write("id\tprecursor_mz\tsize\tidentified_spec_count\tunidentified_spec_count\t"
                  "max_ratio\tmax_il_ratio\tprecursor_mz_range\tsequences\t"
                  "max_sequence\tmax_sequence_count\tmax_sequence_mods\t"
                  "second_max_sequence\tsecond_max_sequence_count\tsecond_max_sequence_mods\tn_input_files")

        if process_synthetic:
            OUT.write("\tsynth_count\tsynth_ratio\tsynth_max_sequence")

        OUT.write("\n")

        # process the file
        parser = clustering_parser.ClusteringParser(input_file)

        for cluster in parser:
            cluster_line = process_cluster(cluster)
            OUT.write(cluster_line)

            # process synthetic peptides
            if process_synthetic:
                synth_line = process_synthetic_peptides(cluster)
                OUT.write("\t" + synth_line)

            OUT.write("\n")

    print("Results written to " + output_file)


if __name__ == "__main__":
    main()
