import os
import sys
import unittest

sys.path.insert(0, os.path.abspath('..'))
import spectra_cluster.ui.cluster_result_comparator as rc


class ClusterAsFeaturesTest(unittest.TestCase):
    def setUp(self):
        self.testfile1 = os.path.join(os.path.dirname(__file__), "testfiles", "kuester_test_short.clustering")
        self.testfile2 = os.path.join(os.path.dirname(__file__), "testfiles", "kuester_test_short_adapted.clustering")

    def testLoadSpectra(self):
        spec1 = rc.load_spectra_to_cluster(self.testfile1)

        self.assertEqual(10, len(spec1))

        # only one cluster should be loaded
        for cluster_id in spec1.values():
            self.assertEqual("7549a1a6-8307-45a5-b351-16182d2daa82", cluster_id)

    def testBuildNetwork(self):
        spec1 = rc.load_spectra_to_cluster(self.testfile1, before_cluster_id="t1")
        spec2 = rc.load_spectra_to_cluster(self.testfile2, before_cluster_id="t2")

        graph = rc.build_network(spec1, spec2, source1="t1", source2="t2")

        self.assertEqual(2, len(graph.nodes))
        self.assertEqual(1, len(graph.edges))

        for e in graph.edges:
            self.assertEqual(10, graph.edges[e]["weight"])

        rc.add_node_properties(graph, self.testfile1, "t1")
        rc.add_node_properties(graph, self.testfile2, "t2")

        # make sure the node properties worked
        for node in graph.nodes:
            self.assertTrue("source" in graph.node[node])
            self.assertEqual(10, graph.node[node]["size"])
            if graph.node[node]["source"] == "source1":
                self.assertEqual("7549a1a6-8307-45a5-b351-16182d2daa82", graph.node[node]["id"])
            else:
                self.assertEqual("8549a1a6-8307-45a5-b351-16182d2daa82", graph.node[node]["id"])