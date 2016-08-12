import unittest
import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import spectra_cluster.analyser.cluster_features as cluster_features
import spectra_cluster.clustering_parser as clustering_parser


class ClusterAsFeaturesTest(unittest.TestCase):
    """
    Test case for the ClusterAsFeatures class
    """
    def setUp(self):
        self.testfile = os.path.abspath('.') + os.path.sep + "test.clustering"

    def testClusterAsFeatures(self):
        parser = clustering_parser.ClusteringParser(self.testfile)

        analyser = cluster_features.ClusterAsFeatures(
            sample_name_extractor=ClusterAsFeaturesTest.pride_project_extractor)

        for cluster in parser:
            analyser.process_cluster(cluster)

        self.assertEqual(838, len(analyser.features))
        self.assertEqual(1, len(analyser.features[0]))
        self.assertEqual(2, analyser.features[0]["PRD000001"])

        self.assertEqual(1, len(analyser.samples))
        self.assertTrue("PRD000001" in analyser.samples)

    @staticmethod
    def pride_project_extractor(spectrum):
        filename = spectrum.get_filename()
        start_index = filename.rfind("/")
        end_index = filename.find(".")

        return filename[start_index + 1:end_index]