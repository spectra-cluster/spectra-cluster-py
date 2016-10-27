#################
Clustering Parser
#################

Usage
=====

The clustering parser provides a simple python iterator::

   import spectra_cluster.clustering_parser as clustering_parser

   clustering_file = "test.clustering"
   parser = clustering_parser.ClusteringParser(clustering_file)

   for cluster in parser:
        # do something with the cluster
        pass

Class Definition
================

.. autoclass:: spectra_cluster.clustering_parser.ClusteringParser
   :members:
