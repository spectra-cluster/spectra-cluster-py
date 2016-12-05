#################
API Documentation
#################

The spectra-cluster-py project provides a set
of APIs that support the development of own
scripts to work with MS/MS clustering results in the
.clustering format.

Core
====

The core of the API is made up of
a set of :doc:`common objects <objects>` that represent
the clustering results and the :doc:`clustering_parser`
used to parse the .clustering files.

Analysers
=========

The functionality of all :doc:`end-user tools <../tools/index>`
is represented through the analysers. Each analyser implements
the :doc:`AbstractAnalyser <abstract_analyser>` class.

Exporter
========

The exporter classes implement the :doc:`AbstractAnalyser <abstract_analyser>`
class but convert the clusters' consensus spectra into different file formats.

##############
API Class List
##############


.. toctree::
   :maxdepth: 2

   objects
   clustering_parser
   abstract_analyser
   id_transferer
   cluster_features
   exporter

