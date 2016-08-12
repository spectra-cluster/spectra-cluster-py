import unittest
import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import spectra_cluster.tools.fasta_paraser as fasta_parser


class FastaParserTest(unittest.TestCase):
    """
    TestCase for the FastaParser class
    """
    def setUp(self):
        self.testfile = os.path.abspath(".") + os.sep + "test.fasta"

    def testIterator(self):
        parser = fasta_parser.FastaParser(self.testfile)

        n_entries = 0

        for fasta_entry in parser:
            n_entries += 1

            if n_entries == 1:
                self.assertEqual(">tr|K7EKG6|K7EKG6_HUMAN Mitochondrial import receptor subunit TOM40 homolog (Fragment) OS=Homo sapiens GN=TOMM40 PE=4 SV=6",
                                 fasta_entry.header_line)
                self.assertEqual("GNVLAASSPPAGPPPPPAPALVGLPPPPPSPPGFT", fasta_entry.sequence)
                self.assertEqual("K7EKG6", fasta_entry.getAccession())

            if n_entries == 198:
                self.assertEqual("ALBU_BOVIN", fasta_entry.getAccession())
            if n_entries == 199:
                self.assertEqual("KKA1_ECOLX", fasta_entry.getAccession())
            if n_entries == 200:
                self.assertEqual("AT1G44020.1", fasta_entry.getAccession())

        self.assertEqual(200, n_entries)