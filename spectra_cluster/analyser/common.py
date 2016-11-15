"""Holds common objects used throughout all analysers"""

import sys


class AbstractAnalyser:
    """Base class for all analysers.

    This class mainly provides helper functions for
    filtering clusters based on size, ratio, identified
    spectra, etc.

    Additionally, every child class must implement the
    **process_cluster(Cluster)** function.

    Its init functions sets the following member variables
    so that all clusters are accepted.

    :ivar min_size: Minimal size of the cluster
    :ivar max_size: Maximum size of the cluster
    :ivar min_ratio: Minimum I/L agnostic ratio
    :ivar max_ratio: Maximum I/L agnostic ratio
    :ivar min_identified_spectra: Minimal identified spectra
    :ivar max_identified_spectra: Maximum identified spectra
    :ivar min_unidentified_spectra: Minimal unidentified spectra
    :ivar max_unidentified_spectra: Maximum unidentified spectra
    """
    def __init__(self):
        """Initialises the default parameters to filter clusters.

        """
        self.min_size = 0
        self.max_size = sys.maxsize
        self.min_ratio = 0
        self.max_ratio = 1
        self.min_identified_spectra = 0
        self.max_identified_spectra = sys.maxsize
        self.min_unidentified_spectra = 0
        self.max_unidentified_spectra = sys.maxsize

    def _ignore_cluster(self, cluster):
        """Tests whether the passed cluster should be ignored

        This function assess the class' instance variables
        (min size, max size, etc.) to test whether the
        passed cluster should be ignored for further
        processing.

        :param: cluster Cluster to test.
        :return: Boolean indicating whether the cluster should be ignored.
        """
        if cluster.n_spectra < self.min_size:
            return True
        if cluster.n_spectra > self.max_size:
            return True
        if cluster.max_il_ratio is None and (self.min_ratio > 0 or self.max_ratio < 1):
            return True
        if cluster.max_il_ratio is not None and cluster.max_il_ratio < self.min_ratio:
            return True
        if cluster.max_il_ratio is not None and cluster.max_il_ratio > self.max_ratio:
            return True
        if cluster.identified_spectra < self.min_identified_spectra:
            return True
        if cluster.identified_spectra > self.max_identified_spectra:
            return True
        if cluster.unidentified_spectra < self.min_unidentified_spectra:
            return True
        if cluster.unidentified_spectra > self.max_unidentified_spectra:
            return True

        return False

    def process_cluster(self, cluster):
        """
        Processes the defined cluster.

        :param cluster: Cluster to process.
        """
        raise NotImplementedError
