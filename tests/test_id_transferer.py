import unittest
import os
from .context import id_transferer
from .context import clustering_parser


class IdTransfererTest(unittest.TestCase):
    """
    Test case for the IdTransfererTest class
    """
    def setUp(self):
        self.testfile = os.path.abspath('.') + os.path.sep + "test.clustering"

    def test_identifiy_all(self):
        parser = clustering_parser.ClusteringParser(self.testfile)
        analyser = id_transferer.IdTransferer(True, True)

        for cluster in parser:
            analyser.process_cluster(cluster)

        self.assertEqual(3151, len(analyser.identification_references))
        self.assertEqual(1, len(analyser.identification_references[11].psms))
        self.assertEqual("HQGVMVGMGQK", analyser.identification_references[11].psms[0].sequence)
        self.assertEqual("/home/jg/Projects/ebi-pride/pride-cluster-2/chimeric-spectra-generator/src/test/resources/PRD000001.st.id.mgf",
                         analyser.identification_references[11].filename)

        ref4 = analyser.identification_references[4]
        self.assertEqual(1, len(ref4.psms))
        self.assertEqual("MEGIGLK", ref4.psms[0].sequence)
