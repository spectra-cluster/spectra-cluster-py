import unittest
import sys
import os

sys.path.insert(0, os.path.abspath('..'))
import spectra_cluster.ui.mgf_search_result_annotator as mgf_search_result_annotator


class MsgfpIdentImportTest(unittest.TestCase):
    """
    Test the MSAmanda import function
    """

    def setUp(self):
        self.ident_file = os.path.join(os.path.dirname(__file__), "testfiles", "msgfplus_test.mzid")
        self.mgf_file = os.path.join(os.path.dirname(__file__), "testfiles", "test.mgf")

    def test_import_testfile(self):
        results = mgf_search_result_annotator.parse_msgfplus_mzident(self.ident_file, 0.05)

        self.assertEqual(len(results), 5, "Incorrect number of PSMs loaded")

        for result in results:
            with self.subTest(result=result):
                self.assertTrue(result.title is not None, msg="Missing title in extracted PSM")


if __name__ == "__main__":
    unittest.main()
