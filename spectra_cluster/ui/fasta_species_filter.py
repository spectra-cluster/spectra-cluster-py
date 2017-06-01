"""fasta_species_filter

This tool only retains proteins from one of the species defined in the passed text file

Usage:
  fasta_species_filter.py --input=<original.fasta> --output=<unique.fasta> --species=<species.txt>

Options:
  -i, --input=<original.fasta>         Path to the FASTA file to process.
  -o, --output=<unique.fasta>          Path to the newly created unique FASTA file
  -s, --species=<species.txt>          Path to a text file containing species mnemonics
                                       (one per line), such as "HUMAN", "BOVIN", etc.
  -h, --help                           Print this help message.
  -v, --version                        Print the current version.
"""

import sys
import os
from docopt import docopt
import hashlib

# make the spectra_cluster packages available
sys.path.insert(0, os.path.abspath('..') + os.path.sep + "..")
from spectra_cluster.tools import fasta_paraser


def main():
    """
    Primary entry function for the CLI.
    :return:
    """
    arguments = docopt(__doc__, version='unique_fasta_extracter 1.0 BETA')

    # make sure the input file exists
    if not os.path.isfile(arguments['--input']):
        print("Error: Cannot find input file '" + arguments["--input"] + "'")
        sys.exit(1)

    if not os.path.isfile(arguments["--species"]):
        print("Error: Cannot find species file '" + arguments["--species"] + "'")
        sys.exit(1)

    # make sure the output file does not exist
    if os.path.isfile(arguments["--output"]):
        print("Error: Output file exists '" + arguments["--output"] + "'")
        sys.exit(1)

    # load the species
    species = list()
    with open(arguments["--species"], "r") as species_in:
        for line in species_in:
            species.append("_" + line.strip().upper())

    parser = fasta_paraser.FastaParser(arguments["--input"])
    total = 0
    retained = 0

    with open(arguments["--output"], "w") as OUT:
        for fasta_entry in parser:
            total += 1

            include_entry = False

            for species_definition in species:
                if species_definition in fasta_entry.header_line:
                    include_entry = True
                    break

            if not include_entry:
                continue

            # write out the entry
            OUT.write(fasta_entry.header_line + "\n")
            OUT.write(fasta_entry.sequence + "\n")

            retained += 1

    print("Wrote " + str(retained) + "/" + str(total) + " sequences to " + arguments["--output"])


if __name__ == "__main__":
    main()
