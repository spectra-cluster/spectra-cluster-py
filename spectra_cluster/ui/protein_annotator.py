"""protein_annotator

This tool adds a protein accession column to a text file based on the present peptides.

The column containing the peptide string must be defined. Then, all proteins that match the given peptide are
written to a specified column.

Usage:
  protein_annotator.py --input=<input.tsv> --output=<extended_file.tsv> --fasta=<fasta_file.fasta>
                       [--peptide_column=<column_name>] [--protein_column=<column_name>]
                       [--protein_separator=<separator>] [--column_separator=<separator>]
                       [--ignore_il]
  protein_annotator.py (--help | --version)

Options:
  -i, --input=<input.tsv>               Path to the input file.
  -o, --output=<extended_file.tsv>      Path to the output file that should be created.
  -f, --fasta=<fasta_file.fasta>        Fasta file to match the peptides to.
  --peptide_column=<column_name>        Column name of the peptide column [default: sequence]
  --protein_column=<column_name>        Column name of the newly added protein column [default: protein]
  --protein_separator=<separator>       Separator to separate multiple protein entries [default: ;]
  --column_separator=<separator>        Separator to separate columns in the file [default: TAB]
  --ignore_il                           If set I/L are treated as synonymous.
  -h, --help                            Print this help message.
  -v, --version                         Print the current version.
"""

import sys
import os
import csv
from docopt import docopt
import re

# make the spectra_cluster packages available
sys.path.insert(0, os.path.abspath('..') + os.path.sep + "..")

from spectra_cluster.tools import fasta_paraser


def extract_separator(user_separator):
    """
    Parses the user defined separator and returns the matching character.

    :param user_separator: The user defined separator.
    :return: The parsed character
    """
    if user_separator == "TAB":
        return "\t"

    return user_separator


def load_peptides(input_file, peptide_column, column_separator):
    """
    Parses the input file and extracts all peptides occuring within the file. Peptide strings
    are cleaned (only valid characters retained) and returned as a set.

    :param input_file: The file to parse
    :param peptide_column: The column header to extract the peptides from.
    :param column_separator: The separator used for the columns
    :return: A set of strings representing the peptides.
    """
    with open(input_file, "r") as input_stream:
        csv_reader = csv.DictReader(input_stream, delimiter=column_separator)

        peptides = set()

        for row in csv_reader:
            if peptide_column not in row:
                raise Exception("Specified peptide column '" + peptide_column + "' not found in input file.")

            sequence = row[peptide_column]
            clean_sequence = re.sub("[^A-Z]", "", sequence)
            peptides.add(clean_sequence)

        return peptides


def map_peptides_to_proteins(peptides, fasta_filename, ignore_il=False):
    """
    Maps the peptides to the proteins in the passed FASTA file.

    :param peptides: A iterable containing the pepitde strings.
    :param fasta_filename: Filename of the FASTA file to parse.
    :param ignore_il: If set to True I/L are treated as the same AA.
    :return: A dict with the peptide as key and the protein accessions as list.
    """
    parser = fasta_paraser.FastaParser(fasta_filename)
    peptide_protein_map = dict()

    for fasta_entry in parser:
        if ignore_il:
            fasta_sequence = fasta_entry.sequence.replace("I", "L")
        else:
            fasta_sequence = fasta_entry.sequence

        # check whether any of the known sequences is present in this protein
        for sequence in peptides:
            # replace all I with L if defined
            if ignore_il:
                comparison_sequence = sequence.replace("I", "L")
            else:
                comparison_sequence = sequence

            if comparison_sequence in fasta_sequence:
                # use the original peptide sequence as a key
                if sequence not in peptide_protein_map:
                    peptide_protein_map[sequence] = list()

                peptide_protein_map[sequence].append(fasta_entry.getAccession())

    return peptide_protein_map


def write_extended_file(input_filename, output_filename, peptides_to_protein, column_separator, protein_separator,
                        peptide_column, protein_column):
    """
    Creates the new file which only is a copy of the current file with the protein column added to the end

    :param input_filename: The input filename path.
    :param output_filename: The output filename path.
    :param peptides_to_protein: A dict containing the peptide string as key and a list of protein acccessions as values.
    :param column_separator: Separator used for the columns.
    :param protein_separator: Separator used if multiple protein accessions are found.
    :param peptide_column: Name of the peptide column.
    :param protein_column: New name of the protein column.
    :return:
    """
    # get the index of the peptide column
    peptide_index = -1
    first_line = True
    with open(output_filename, "w") as output_file:
        with open(input_filename, "r") as input_file:
            for line in input_file:
                # remove the trailing line separator
                line = line.rstrip("\n")

                if first_line:
                    new_line = line + column_separator + protein_column + "\n"
                    output_file.write(new_line)
                    first_line = False

                    # get the index of the peptide column
                    fields = line.split(column_separator)
                    for i in range(0, len(fields)):
                        if fields[i] == peptide_column:
                            peptide_index = i
                            break
                    if peptide_index < 0:
                        raise Exception("Failed to find peptide column '" + peptide_column + "' in input file.")
                else:
                    # get the peptide string
                    fields = line.split(column_separator)
                    sequence = fields[peptide_index]
                    clean_sequence = re.sub("[^A-Z]", "", sequence)

                    # write the original line
                    output_file.write(line + column_separator)

                    if clean_sequence in peptides_to_protein and peptides_to_protein[clean_sequence] is not None:
                        output_file.write(protein_separator.join(peptides_to_protein[clean_sequence]))

                    output_file.write("\n")


class ProteinMappings:
    """
    Simple class to represent all proteins a peptide maps to
    """
    def __init__(self, proteins):
        """
        Creates a new PeptideToproteinMapping object
        :param proteins: A list of proteins the pepitde maps to.
        :return:
        """
        self.proteins = proteins
        self.n_proteins = len(self.proteins)


class ProteinGroup:
    """
    Simple representation of a ProteinGroup
    """
    def __init__(self, label, accessions):
        """
        Creates a new protein group object

        :param label: The label of the protein group.
        :param accessions: The accessions of all member proteins
        :return:
        """
        self.label = label
        self.accessions = tuple(accessions)

    def __hash__(self):
        return hash(self.label)


def do_protein_inference(peptides_to_proteins, protein_separator=";"):
    """
    Warning: This function currently does **NOT** work correctly!

    Creates the smallest set of protein (groups) required to explain all peptides. Ambiguous peptides
    are not mapped to any group.

    :param peptides_to_proteins: A dict with the peptide sequence as key and all mapping proteins as a list (value).
    :param protein_separator: The separator to use when creating the label for protein groups.
    :return: A dict with the peptide sequence as key and the protein / protein group accession as value (single entry
             in a list for compatibility reasons).
    """
    raise ArithmeticError()

    # create list of PeptideToProteinMappings
    mappings = [ProteinMappings(proteins) for proteins in peptides_to_proteins.values()]

    # extract the smallest set of protein groups necessary to explain all peptides

    # start with the peptides mapping to the smallest number of proteins (ie. partially unique ones)
    mappings.sort(key=lambda x: x.n_proteins, reverse=False)

    retained_proteins = set()
    retained_protein_groups = list()

    for mapping in mappings:
        # if one of the accessions is already retained, simply ignore the rest
        if any(True for accession in mapping.proteins if accession in retained_proteins):
            continue

        # create a new protein group
        retained_proteins |= set(mapping.proteins)
        retained_protein_groups.append(ProteinGroup(protein_separator.join(mapping.proteins), mapping.proteins))

    # create a map of proteins to protein groups
    protein_to_protein_group = dict()

    for protein_group in retained_protein_groups:
        for protein in protein_group.accessions:
            protein_to_protein_group[protein] = protein_group

    # remap the peptides
    new_peptide_mappings = dict()

    for sequence in peptides_to_proteins:
        all_proteins = peptides_to_proteins[sequence]

        # get all matching protein groups
        matching_protein_groups = set()

        for protein in all_proteins:
            # if the protein is a subset protein, it was ignored
            if protein not in protein_to_protein_group:
                continue

            matching_protein_groups.add(protein_to_protein_group[protein])

        # if the peptide matches multiple protein groups, ignore it
        if len(matching_protein_groups) > 1:
            new_peptide_mappings[sequence] = None
        elif len(matching_protein_groups) == 1:
            new_peptide_mappings[sequence] = [list(matching_protein_groups)[0].label]
        else:
            # missing mappings are also represented by None
            new_peptide_mappings[sequence] = None

    return new_peptide_mappings


def main():
    """
    Primary entry function for the CLI.
    """
    arguments = docopt(__doc__, version='protein_annotator 1.0 BETA')

    # make sure the input files exist
    if not os.path.isfile(arguments['--input']):
        print("Error: Cannot find input file '" + arguments["--input"] + "'")
        sys.exit(1)
    if not os.path.isfile(arguments['--fasta']):
        print("Error: Cannot find fasta file '" + arguments["--fasta"] + "'")
        sys.exit(1)

    # make sure the output file does not exist
    if os.path.isfile(arguments["--output"]):
        print("Error: Output file exists '" + arguments["--output"] + "'")
        sys.exit(1)

    # get the separator
    column_separator = extract_separator(arguments["--column_separator"])
    protein_separator = extract_separator(arguments["--protein_separator"])

    peptide_column = arguments["--peptide_column"]
    protein_column = arguments["--protein_column"]

    # load all peptides
    print("Loading peptides from file...", end="")
    peptides = load_peptides(arguments["--input"], peptide_column, column_separator)
    print("Done. (" + str(len(peptides)) + " loaded)")

    # map the proteins
    print("Mapping peptides to proteins...", end="")
    peptides_to_protein = map_peptides_to_proteins(peptides, arguments["--fasta"], arguments["--ignore_il"])
    print("Done.")

    # if arguments["--protein_inference"]:
    #    peptides_to_protein = do_protein_inference(peptides_to_protein, protein_separator)

    # write the new file
    write_extended_file(arguments["--input"], arguments["--output"], peptides_to_protein, column_separator,
                        protein_separator, peptide_column, protein_column)

if __name__ == "__main__":
    main()
