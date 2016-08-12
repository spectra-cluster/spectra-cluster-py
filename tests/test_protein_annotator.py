import unittest
import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import spectra_cluster.ui.protein_annotator as protein_annotator


class FastaParserTest(unittest.TestCase):
    """
    TestCase for the ProteinParser
    """
    def setUp(self):
        self.testfile = os.path.abspath(".") + os.sep + "test.fasta"

    def testMapPeptides(self):
        peptides = set()
        peptides.add("KWVTFISLLLL")
        peptides.add("AGGE")
        peptides.add("LLA")
        peptides.add("A")

        protein_map = protein_annotator.map_peptides_to_proteins(peptides, self.testfile)

        self.assertEqual(len(peptides), len(protein_map))

        for sequence in protein_map.keys():
            org_len = len(protein_map[sequence])
            unique_len = len(set(protein_map[sequence]))

            self.assertEqual(unique_len, org_len)

    def testProteinInference(self):
        peptides = set()
        peptides.add("GLL")
        peptides_to_proteins = protein_annotator.map_peptides_to_proteins(peptides, self.testfile)
        peptide_map = protein_annotator.do_protein_inference(peptides_to_proteins)

        self.assertEqual(1, len(peptide_map))
        self.assertTrue("GLL" in peptide_map)
        self.assertEqual(24, peptide_map["GLL"][0].count(";"))

        # add another unique peptide
        peptides.add("EEWQCLDTAQRNLYKNV")
        peptides_to_proteins = protein_annotator.map_peptides_to_proteins(peptides, self.testfile)
        peptide_map = protein_annotator.do_protein_inference(peptides_to_proteins)

        self.assertEqual(2, len(peptide_map))
        self.assertEqual("M0QXM7", peptide_map["EEWQCLDTAQRNLYKNV"][0])
        self.assertEqual("M0QXM7", peptide_map["GLL"][0])

    def testIgnoreIL(self):
        peptides = set()
        peptides.add("GLL")

        peptides_to_proteins = protein_annotator.map_peptides_to_proteins(peptides, self.testfile, ignore_il=False)
        self.assertEqual(25, len(peptides_to_proteins["GLL"]))

        peptides_to_proteins = protein_annotator.map_peptides_to_proteins(peptides, self.testfile, ignore_il=True)
        self.assertEqual(46, len(peptides_to_proteins["GLL"]))

