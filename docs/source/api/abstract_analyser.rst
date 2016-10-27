#################
Abstract Analyser
#################

Usage
=====

This class provides some common features often required in analysers. Currently, the main feature is
the filtering of clusters based on size, ratio, identified spectra, etc.

If you develop an analyser simply create a child of this class::

    from spectra_cluster.analyser.common import AbstractAnalyser

    class MyAnalyser(AbstractAnalyser):
        """
            This Analyser counts the total number of
            clusters
        """
        def __init__(self):
            # call the AbstractAnalyser's init class
            # to set all filtering variables to their
            # default values.
            super().__init__()

            # initialise you're class' instance
            # variables            self.number_of_clusters = 0

        def process_cluster(self, cluster):
            # This function must be implemented by all
            # analysers

            # first, check whether this cluster should
            # be processed at all
            if self._ignore_cluster(cluster):
                return

            # count the cluster
            self.number_of_clusters += 1

Class Definition
================

.. autoclass:: spectra_cluster.analyser.common.AbstractAnalyser
   :members:
   :private-members: