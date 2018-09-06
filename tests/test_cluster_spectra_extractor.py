import unittest
import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import spectra_cluster.ui.cluster_spectra_extractor as cluster_spectra_extractor
import pickle


class ClusterSpectraExtractorTest(unittest.TestCase):
    def setUp(self):
        self.testfile = os.path.join(os.path.dirname(__file__), "testfiles", "msamanda_test_output.mgf")

    def test_mgf_index(self):
        cluster_spectra_extractor.build_mgf_indices([self.testfile])

        index_file = self.testfile + ".pyindex"

        self.assertTrue(os.path.isfile(index_file), msg="MGF index file was not created")

        # load the index file
        with open(index_file, "rb") as reader:
            index = pickle.load(reader)

        self.assertEqual(51, len(index))
        self.assertEqual(0, index[0])

        with open(self.testfile, "r") as mgf_reader:
            for spec_index in index:
                mgf_reader.seek(spec_index)
                line = mgf_reader.readline()

                self.assertEqual("BEGIN IONS", line.strip())

        os.remove(index_file)
        # TODO: test append_spectra_to_file with and without index
