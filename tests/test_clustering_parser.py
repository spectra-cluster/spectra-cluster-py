import unittest
import os
from .context import clustering_parser


class ClusteringParserTest(unittest.TestCase):
    """
    Test case for the ClusteringParser class
    """
    def setUp(self):
        self.testfile = os.path.abspath('.') + os.path.sep + "test.clustering"
        self.ptm_string = "1-MOD:1234,2-MOD:00043"
        self.spec_line = "SPEC\t#file=/home/jg/Projects/ebi-pride/pride-cluster-2/chimeric-spectra-generator/src/test/" \
                         "resources/PRD000001.st.id.mgf#id=index=1464#title=id=PRD000001;PRIDE_Exp_Complete_Ac_1644.xml;" \
                         "spectrum=5071,splib_sequence=MEGIGLK,score=0.362,peptideR2=,scoreR2=\ttrue\tMEGIGLK\t382.149\t2" \
                         "\t\t\t0.0"

    def test_parse_ptms(self):
        ptms = clustering_parser.ClusteringParser._parse_ptms(self.ptm_string)

        self.assertEqual(2, len(ptms))
        self.assertEqual(1, ptms[0].position)
        self.assertEqual("MOD:1234", ptms[0].accession)
        self.assertEqual(2, ptms[1].position)
        self.assertEqual("MOD:00043", ptms[1].accession)

    def test_parse_spec_line(self):
        psm = clustering_parser.ClusteringParser._parse_spec_line(self.spec_line)

        self.assertEqual("MEGIGLK", psm.sequences[0])
        self.assertEqual(2, psm.charge)
        self.assertEqual(382.149, psm.precursor_mz)

        self.assertEqual("/home/jg/Projects/ebi-pride/pride-cluster-2/chimeric-spectra-generator/src/test/resources/"
                         "PRD000001.st.id.mgf", psm.get_filename())
        self.assertEqual("index=1464", psm.get_id())
        self.assertEqual("id=PRD000001;PRIDE_Exp_Complete_Ac_1644.xml;spectrum=5071,splib_sequence=MEGIGLK,"
                         "score=0.362,peptideR2=,scoreR2=", psm.get_title())

        self.assertEqual("MEGIGLK", psm.get_clean_sequences()[0])
        psm.sequences = ["AsdEs(+32)DE"]
        self.assertEqual("ASDESDE", psm.get_clean_sequences()[0])

    def test_parse_clustering_file(self):
        parser = clustering_parser.ClusteringParser(self.testfile)

        n_clusters = 0

        for cluster in parser:
            n_clusters += 1

            if n_clusters == 1:
                self.assertEqual("1cc813a1-4e75-4c1d-99aa-752312fbe554", cluster.id)
                self.assertEqual(359.155, cluster.precursor_mz)
                self.assertEqual(2, len(cluster.spectra))
                self.assertEqual(1, len(cluster.spectra[0].sequences))
                self.assertEqual("RPHFFFPK", cluster.spectra[0].sequences[0])

        self.assertEqual(838, n_clusters)
