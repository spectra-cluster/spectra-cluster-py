#################################################
spectra-cluster-py - Analysing clustering results
#################################################

Welcome
=======

The spectra-cluster-py project is a collection of tools
and APIs that help analysing and working with MS/MS spectrum
clustering results in the **.clustering format**.

    * **Download** the
      `latest release <https://github.com/spectra-cluster-py/releases>`_
    * **Report problems** using the project's
      `issue tracker <https://github.com/spectra-cluster/spectra-cluster-py/issues>`_
    * **Fork** the project on GitHub
      https://github.com/spectra-cluster/spectra-cluster-py

The .clustering format
======================

The .clustering format is currently used by the
`spectra-cluster`_ algorithm (and API) and the
output format of the derived tools:

    * `spectra-cluster-cli`_: Command line version
      of the `spectra-cluster`_ algorithm.
    * `spectra-cluster-gui`_: Graphical user interface
      for running the `spectra-cluster-cli`_ tool.
    * `spectra-cluster-hadoop`_: A Hadoop implementation
      of the `spectra-cluster`_ algorithm.

An up-to-date documentation of the .clustering format can
be found at the Java API `clustering-file-reader`_ project
page.

.. _spectra-cluster: https://github.com/spectra-cluster/spectra-cluster
.. _spectra-cluster-cli: https://github.com/spectra-cluster/spectra-cluster-cli
.. _spectra-cluster-gui: https://github.com/spectra-cluster/spectra-cluster-gui
.. _spectra-cluster-hadoop: https://github.com/spectra-cluster/spectra-cluster-hadoop
.. _clustering-file-reader: https://github.com/clustering-file-reader

Tools
=====

The spectra-cluster-py project contains a set of end-user
ready tools to analyse MS/MS clustering results in the
.clustering format.

The :doc:`id_transferer_cli <tools/id_transferer_cli>`
transfers identification data to
unidentified spectra. This can be used to improve
the accuracy of **label-free** quantitation data.

The :doc:`cluster_features_cli <tools/cluster_features_cli>`
creates a matrix with the input MGF file names as column
headers and the clusters as rows. Each cell contains the number
of spectra per file and cluster. This can be used, for example, to
run a principal component analysis of the input files based
on the clusters.

The :doc:`protein_annotator <tools/protein_annotator>` can map
peptides in a text file to proteins from a fasta file. Additionally,
basic protein inference can be performed.

The complete list of tools can be found
:doc:`here <tools/index>`.

In the python package the source code of these tools is
at ``spectra_cluster.ui``.

Python API
==========

This collection of classes is intended to help you develop your
own scripts to analyse MS/MS clustering results in the
**.clustering** format.

You can find the complete API documentation
:doc:`here <api/index>`.

Documentation
=============

* :doc:`contents` for a complete overview
* :ref:`Search page <search>` search the documentation
* :ref:`genindex`
* :ref:`modindex`


