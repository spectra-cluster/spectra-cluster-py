import unittest
import sys
import os
from pyteomics import mzid

sys.path.insert(0, os.path.abspath('..'))
import spectra_cluster.ui.mgf_search_result_annotator as mgf_search_result_annotator


class MzIdParserTest(unittest.TestCase):
    """
    Test mzid parser functions
    """
    def setUp(self):
        self.testfile = os.path.join(os.path.dirname(__file__), "testfiles", "test.mzid")

    def testPyteomicsParsing(self):
        return
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

    def testParseMzIdentMl(self):
        return
        score_field =mgf_search_result_annotator.get_scorefield_mzident(filename=self.testfile)

        self.assertEqual("Scaffold:Peptide Probability", score_field)
        resutls = mgf_search_result_annotator.parser_mzident(filename=self.testfile,
                                                             score_field=score_field,
                                                             fdr=2)
        self.assertEqual(6573, len(resutls))

    def testGetPeakFile(self):
        return
        peakfile =mgf_search_result_annotator.get_source_peak_file_mzident(filename=self.testfile)
        self.assertEqual("R2_F53 (R2_F53_F094867).mzid_R2_F53_(R2_F53_F094867).MGF", peakfile)


    def testGetScfieldPeakfile(self):
        return
        score_field, peakfile = mgf_search_result_annotator.get_scfield_peakfile(filename=self.testfile)
        self.assertEqual("Scaffold:Peptide Probability", score_field)
        self.assertEqual("R2_F53 (R2_F53_F094867).mzid_R2_F53_(R2_F53_F094867).MGF", peakfile)


if __name__ == "__main__":
    unittest.main()
