import unittest
import sys
import os
from pyteomics import mzid

sys.path.insert(0, os.path.abspath('..'))
import spectra_cluster.ui.mgf_search_result_annotator as mgf_search_result_annotator


class ScaffoldImportTest(unittest.TestCase):
    """
    Test scaffold mzid import functions
    """
    def setUp(self):
        self.testfile = os.path.sep.join([os.path.abspath("."), "testfiles", "test.mzid"])

    def testPyteomicsParsing(self):
        reader = mzid.read(self.testfile)

        n_decoy = 0
        n_target = 0
        for spec_ref in reader:
            for spec_id in spec_ref["SpectrumIdentificationItem"]:
                if spec_id["passThreshold"]:
                    n_target += 1
                else:
                    n_decoy += 1

            # spectrumID
            # spectrum title
            # SpectrumIdentificationItem (list)
            # - dict
            # - Mascot:score
            # - Mascot:identify threshold
            # - Scaffold:Peptide Probability
            # - chargeState
            # - experimentalMassToCharge
            # - passThreshold (bool)
            # - peptide_ref (str = sequence)
            # - rank (int)

        self.assertEqual(6573, n_target)
        self.assertEqual(14095, n_decoy)

    def testParseScaffold(self):
        results = mgf_search_result_annotator.parse_scaffold(self.testfile, 2)

        self.assertEqual(6573, len(results))

    def testParseMzIdentMl(self):
        resutls = mgf_search_result_annotator.parser_mzident(filename=self.testfile,
                                                             score_field="Scaffold:Peptide Probability",
                                                             fdr=2)
        self.assertEqual(6573, len(resutls))
