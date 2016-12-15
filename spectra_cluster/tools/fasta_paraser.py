"""
fasta_parser provides the FastaParser class. It allows one to iterate
over all entries within a FASTA file. Fasta entries are available through
the FastaEntry class.
"""


class FastaParser:
    """
    Class to iterate over FASTA files
    """
    def __init__(self, fasta_filename):
        self.fasta_filename = fasta_filename

    def _get_iterator(self):
        """
        Creates an iterator over all entries within the FASTA field
        using yield

        :yield: FastaEntry objects
        """
        with open(self.fasta_filename, "r") as in_file:
            current_header = None
            current_sequence = list()

            for line in in_file:
                line = line.strip()

                if len(line) == 0:
                    continue

                if line[0] == ">":
                    if current_header is not None:
                        entry = FastaEntry(current_header, "".join(current_sequence))
                        yield entry

                    current_header = line
                    current_sequence = list()
                else:
                    current_sequence.append(line)

            # return the final entry
            if current_header is not None:
                entry = FastaEntry(current_header, "".join(current_sequence))
                yield entry

    def __iter__(self):
        """
        Returns an iterator over all fasta entries.

        :return:
        """
        return self._get_iterator()

class FastaEntry:
    """
    Represents a fasta entry
    """
    def __init__(self, header_line, sequence):
        self.header_line = header_line
        self.sequence = sequence

    def containsSubsequence(self, sequence):
        """
        Tests whether the passed sequence is part of the FASTA sequence

        :param sequence: The sub-sequence (ie peptide) to test.
        :return: boolean whether the sequence in a subsequence of the entry.
        """
        return sequence in self.sequence

    def getAccession(self, format="detect"):
        """
        Retrieves the accession for the current FASTA entry.

        :param format: The header line format to expect. If set to "detect" the function tries to guess the format.
                       The following formats are supported "detect", "uniprot", "first_word"
        :return:
        """
        if format == "detect":
            format = FastaEntry.detect_header_format(self.header_line)

        if format == "uniprot":
            start = 4
            end = self.header_line.find("|", start)

        if format == "first_word":
            start = 1
            end_blank = self.header_line.find(" ")
            end_pipe = self.header_line.find("|")

            if end_blank < end_pipe:
                end = end_blank
            else:
                end = end_pipe

            if end < start:
                end = len(self.header_line)

        return self.header_line[start:end]

    @staticmethod
    def detect_header_format(header_line):
        """
        Detects the header format and returns the possible format name.

        :param header_line: The header line to analyse
        :return: The format name ("uniprot", "first_word")
        """
        if header_line[0:4] == ">sp|" or header_line[0:4] == ">tr|" or header_line[0:4] == ">up|":
            return "uniprot"

        return "first_word"
