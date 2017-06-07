import unittest
import os
import sys
import tempfile
sys.path.insert(0, os.path.abspath('..'))
import spectra_cluster.analyser.cluster_features as cluster_features
import spectra_cluster.clustering_parser as clustering_parser


class ClusterAsFeaturesTest(unittest.TestCase):
    """
    Test case for the ClusterAsFeatures class
    """
    def setUp(self):
        self.testfile = os.path.join(os.path.dirname(__file__), "test.clustering")

    def testClusterAsFeatures(self):
        result_file = tempfile.TemporaryFile(mode="w+")

        parser = clustering_parser.ClusteringParser(self.testfile)

        analyser = cluster_features.ClusterAsFeatures(
            result_file=result_file,
            sample_name_extractor=ClusterAsFeaturesTest.pride_project_extractor)

        for cluster in parser:
            analyser.process_cluster(cluster)

        # read the result
        features = list()
        result_file.seek(0)

        for line in result_file:
            fields = line.split("\t")
            spectra_per_cluster = dict()

            for i in range(1, len(fields)):
                sample_id = analyser.sample_ids[i - 1]
                spectra_per_cluster[sample_id] = int(fields[i])

            features.append(spectra_per_cluster)

        self.assertEqual(838, len(features))
        self.assertEqual(1, len(features[0]))
        self.assertEqual(2, features[0]["PRD000001"])

        self.assertEqual(1, len(analyser.sample_ids))
        self.assertTrue("PRD000001" in analyser.sample_ids)

        result_file.close()

    @staticmethod
    def pride_project_extractor(spectrum):
        filename = spectrum.get_filename()
        start_index = filename.rfind("/")
        end_index = filename.find(".")

        return filename[start_index + 1:end_index]

if __name__ == "__main__":
    unittest.main()
