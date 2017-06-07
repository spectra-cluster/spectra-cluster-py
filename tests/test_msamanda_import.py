import unittest
import sys
import os

sys.path.insert(0, os.path.abspath('..'))
import spectra_cluster.ui.mgf_search_result_annotator as mgf_search_result_annotator


class MsAmandaImportTest(unittest.TestCase):
    """
    Test the MSAmanda import function
    """
    def setUp(self):
        self.msamanda_file = os.path.join(os.path.dirname(__file__), "testfiles", "msamanda_test_output.csv")
        self.mgf_file = os.path.join(os.path.dirname(__file__), "testfiles", "msamanda_test_output.mgf")
        self.org_mgf_file = "20151023_QexHF1_RSLC5_Elling_Elling_IMBA_shotgun_TMT_SCX_30per_fr105.mgf"

        # load the MGF titles
        self.mgf_titles = list()

        with open(self.mgf_file, "r") as INPUT_FILE:
            for line in INPUT_FILE:
                line = line.strip()

                if line[0:6] == "TITLE=":
                    self.mgf_titles.append(line[:6])

    def test_import_testfile(self):
        results = mgf_search_result_annotator.parse_msamanda(self.msamanda_file, 0.05, self.org_mgf_file)

        self.assertEqual(len(results), 4, "Incorrect number of PSMs loaded")

        for result in results:
            with self.subTest(result=result):
                self.assertTrue(result.title is not None, msg="Missing title in extracted PSM")

if __name__ == "__main__":
    unittest.main()
