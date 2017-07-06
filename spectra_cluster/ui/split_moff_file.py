"""split_moff_file

Split a MoFF quantification output file based on the MGF input file that the
various PSMs originate from.

Usage:
  split_moff_file.py --input=<moff_result.txt> --out_dir=</my/dir>
  split_moff_file.py (--help | --version)

Options:
  -i, --input=<moff_result.txt>        The path to the moff result file. New files will be generated
                                       by appending "moff.txt" to the MGF filename.
  -o, --out_dir=</my/dir>              Output directory to save the newly created files in.
  -h, --help                           Print this help message.
  -v, --version                        Print the current version.
"""

import sys
import os
from docopt import docopt


def main():
    arguments = docopt(doc=__doc__, version='split_moff_file.py 1.0 ALPHA')

    in_file = arguments["--input"]
    out_dir = arguments["--out_dir"]

    if not os.path.isfile(in_file):
        print("Error: Cannot find input file: " + in_file)
        sys.exit(1)

    if not os.path.isdir(out_dir):
        print("Error: Output directory '" + out_dir + "' does not exist")
        sys.exit(1)

    out_file_handels = dict()
    with open(in_file, "r") as reader:
        header = next(reader)

        # get the column index of the filename column
        filename_index = -1
        header_fields = header.split("\t")
        for i in range(0, len(header_fields)):
            if header_fields[i] == "filename":
                filename_index = i
                break

        if filename_index < 0:
            print("Error: Cannot find 'filename' column.\n")
            sys.exit(1)

        for line in reader:
            fields = line.strip().split("\t")

            filename = os.path.basename(fields[filename_index])

            # open the file if it hasn't been opened yet and write the header line
            if filename not in out_file_handels:
                out_file_handels[filename] = open(os.path.join(out_dir, filename + ".moff.txt"), "w")
                out_file_handels[filename].write(header)

            # write out the line
            out_file_handels[filename].write(line)

    # close all out files
    for filename in out_file_handels:
        out_file_handels[filename].close()

    print(str(len(out_file_handels)) + " new files created.")


if __name__ == "__main__":
    main()
