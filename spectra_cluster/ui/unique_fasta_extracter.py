"""unique_fasta_extracter

This tool simply removes all duplicate protein entries (based on the sequence)
from a given FASTA file.

Usage:
  unique_fasta_extracter.py --input=<original.fasta> --output=<unique.fasta>

Options:
  -i, --input=<original.fasta>         Path to the FASTA file to process.
  -o, --output=<unique.fasta>          Path to the newly created unique FASTA file
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

    # make sure the output file does not exist
    if os.path.isfile(arguments["--output"]):
        print("Error: Output file exists '" + arguments["--output"] + "'")
        sys.exit(1)

    parser = fasta_paraser.FastaParser(arguments["--input"])
    known_sequences = set()
    total = 0
    retained = 0

    with open(arguments["--output"], "w") as OUT:
        for fasta_entry in parser:
            total += 1

            # create the hash
            sequence_hash = hashlib.sha256(fasta_entry.sequence.encode()).hexdigest()

            if sequence_hash in known_sequences:
                continue

            # write out the entry
            OUT.write(fasta_entry.header_line + "\n")
            OUT.write(fasta_entry.sequence + "\n\n")

            # remember the hash
            known_sequences.add(sequence_hash)

            retained += 1

    print("Wrote " + str(retained) + "/" + str(total) + " sequences to " + arguments["--output"])


if __name__ == "__main__":
    main()