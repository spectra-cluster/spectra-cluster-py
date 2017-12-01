"""cluster_spectra_extractor

This tool extracts the spectra of a specified cluster from the original source peak
list files and writes all of a cluster's spectra in a single MGF file.

Usage:
  cluster_spectra_extractor --output_directory=</path/to/results> --clustering_file=<result.clustering>
                            [--add_consensus_spectrum] --peaklist_directory=</path/to/dir>... <cluster_id>...
  cluster_spectra_extractor --build_index <mgf_files>...
  cluster_spectra_extractor (--help | --version)

Options:
  -o, --output_directory=</path/to/results>        Path to the directory where the MGF files will be created. The
                                                   files will have the cluster's id as a name.
  -c, --clustering_file=<result.clustering>        Path to the .clustering result file.
  -p, --peaklist_directory=</path/to/dir>          Path to a directory holding the original MGF files. Multiple
                                                   directories can be specified by specifying this parameter multiple
                                                   times.
  -i, --build_index                                If set, the passed MGF files are indexed.
  --add_consensus_spectrum                         If set, the cluster's consensus spectrum is written to the MGF file
                                                   as the first spectrum.
  -h, --help                                       Displays this help.
  -v, --version                                    Displays the tool's version.
"""

import os
import sys
import pickle
from docopt import docopt
from spectra_cluster import clustering_parser


# This list is used to make sure that additional parameters
# are added in the same order to spectra (if they are available)
possible_mgf_params = ["MIN_COMP", "ADDING_SCORE"]


class SpectrumReference:
    def __init__(self, index, params):
        self.index = index
        self.params = params


class MgfFile:
    """
    Represents an MGF file that can be indexed. This class does not
    represent the actual content of an MGF file but only supports
    random retrieval of spectra
    """
    def __init__(self, filename):
        """
        Create a new MgfFile instance
        :param filename: Path to the MGF file
        """
        self.filename = filename
        self.index_file = filename + ".pyindex"
        self.is_indexed = os.path.isfile(self.index_file)

        if self.is_indexed:
            with open(self.index_file, "rb") as reader:
                self.index = pickle.load(reader)
        else:
            self.index = None

    def get_spectrum_string(self, spec_index):
        """
        Load the spectrum (as a string) defined by the 1-based index.
        :param spec_index: 1-based index of the spectrum to load
        :return: A list of strings each representing one line of the spectrum.
        """
        if spec_index < 1:
            raise Exception("MGF spectrum indices are 1-based")

        # get the offset
        if self.index:
            if spec_index > len(self.index):
                raise Exception("MGF file only contains " + str(len(self.index)) + " spectra. Cannot "
                                                                                   "load spectrum " +
                                str(spec_index))
            offset = self.index[spec_index - 1]
        else:
            with open(self.filename, "r") as reader:
                cur_offset = reader.tell()
                cur_spec = 0
                cur_line = reader.readline()

                while cur_line:
                    if cur_line[:10] == "BEGIN IONS":
                        cur_spec += 1

                    cur_offset = reader.tell()
                    cur_line = reader.readline()

                    if cur_spec == spec_index:
                        offset = cur_offset
                        break

        if not offset:
            raise Exception("Failed to find spectrum " + str(spec_index) + " in " + os.path.basename(self.filename))

        with open(self.filename, "r") as reader:
            reader.seek(offset)
            spec_lines = list()
            line = reader.readline()

            while line:
                spec_lines.append(line)
                if line[:8] == "END IONS":
                    return spec_lines

                line = reader.readline()


def extract_clusters(cluster_ids, clustering_file):
    parser = clustering_parser.ClusteringParser(clustering_file)
    clusters = list()

    for cluster in parser:
        if cluster.id in cluster_ids:
            clusters.append(cluster)
        if len(clusters) == len(cluster_ids):
            break

    return clusters


def get_spectra_per_file(cluster):
    spec_per_file = dict()

    for spectrum in cluster.get_spectra():
        mgf_file = os.path.basename(spectrum.get_filename())

        spec_ref = SpectrumReference(spectrum.get_id(), spectrum.properties)

        if mgf_file not in spec_per_file:
            spec_per_file[mgf_file] = [spec_ref]
        else:
            spec_per_file[mgf_file].append(spec_ref)

    return spec_per_file


def append_spectra_to_file(in_file, spec_refs, out_file):
    with open(out_file, "a") as writer:
        mgf_file = MgfFile(in_file)
        for spec_ref in spec_refs:
            spec_lines = mgf_file.get_spectrum_string(int(spec_ref.index[6:]))

            last_param_line = 0
            for line in spec_lines:
                if line[0].isnumeric():
                    break
                last_param_line += 1

            # write the params first
            writer.write("".join(spec_lines[0:last_param_line]))

            # write any additional params
            for param_name in possible_mgf_params:
                if param_name in spec_ref.params:
                    writer.write(param_name + "=" + spec_ref.params[param_name] + "\n")

            writer.write("".join(spec_lines[last_param_line:]) + "\n")


def build_mgf_indices(mgf_files):
    for mgf_file in mgf_files:
        print("Indexing " + mgf_file + "...")

        # create the index
        index = list()

        with open(mgf_file, "r") as reader:
            cur_offset = reader.tell()
            cur_line = reader.readline()

            while cur_line:
                # simply store the start of a spectrum
                if cur_line[:10] == "BEGIN IONS":
                    index.append(cur_offset)

                cur_offset = reader.tell()
                cur_line = reader.readline()

            # save the index as a pickle object
            with open(mgf_file + ".pyindex", "wb") as writer:
                pickle.dump(index, file=writer)


def write_consensus_spectrum(cluster, mgf_file):
    """
    Writes the cluster's consensus spectrum to the specified (MGF) file
    :param cluster: The cluster object
    :param mgf_file: Path to the result file. This file will be overwritten.
    """
    with open(mgf_file, "w") as writer:
        writer.write("BEGIN IONS\nTITLE=Consensus;cluster_id=" + cluster.id + "\n")
        writer.write("PEPMASS=" + str(cluster.precursor_mz) + "\n")
        writer.write("CHARGE=" + str(cluster.charge) + "+\n")

        for i in range(0, len(cluster.consensus_mz)):
            writer.write(str(cluster.consensus_mz[i]) + " " + str(cluster.consensus_intens[i]) + "\n")

        writer.write("END IONS\n\n")


def main():
    args = docopt(__doc__, version="cluster_spectra_extractor 1.0 BETA")

    if args["--build_index"]:
        build_mgf_indices(args["<mgf_files>"])
        return

    out_dir = args["--output_directory"]
    clustering_file = args["--clustering_file"]
    peak_dirs = args["--peaklist_directory"]
    cluster_ids = args["<cluster_id>"]
    add_consensus = args["--add_consensus_spectrum"]

    if not os.path.isdir(out_dir):
        print("Error: Output directory does not exist")
        sys.exit(1)
    if not os.path.isfile(clustering_file):
        print("Error: Cannot find .clustering file")
        sys.exit(1)
    for peak_dir in peak_dirs:
        if not os.path.isdir(peak_dir):
            print("Error: " + peak_dir + " does not exist")
            sys.exit(1)

    # get the cluster objects for every cluster id
    print("Loading clusters from .clustering file...")
    clusters = extract_clusters(cluster_ids, clustering_file)

    if len(clusters) < len(cluster_ids):
        print("Warning: only " + str(len(clusters)) + " found out of " + str(len(cluster_ids)))
    else:
        print("  " + str(len(clusters)) + " clusters loaded")

    # process the clusters
    for cluster in clusters:
        print("Processing " + cluster.id + "...")

        # get the spectra per mgf file
        spectra_per_file = get_spectra_per_file(cluster)

        output_name = os.path.join(out_dir, cluster.id + ".mgf")

        # if set, add the cluster's consensus spectrum
        if add_consensus:
            write_consensus_spectrum(cluster, output_name)

        for mgf_filename in spectra_per_file.keys():
            # check the directories whether they contain the file
            complete_name = None

            for peak_dir in peak_dirs:
                if os.path.isfile(os.path.join(peak_dir, mgf_filename)):
                    complete_name = os.path.join(peak_dir, mgf_filename)
                    break

            if complete_name is None:
                print("Error: Failed to find " + mgf_filename)
                sys.exit(1)

            # write the spectra
            print("  Extracting " + str(len(spectra_per_file[mgf_filename])) + " spectra from " + mgf_filename + "...")
            append_spectra_to_file(complete_name, spectra_per_file[mgf_filename], output_name)


if __name__ == "__main__":
    main()
