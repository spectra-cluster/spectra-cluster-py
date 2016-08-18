"""The IdTransferer analyser transfers identifications to spectra part of a cluster.

The analysis is run by calling 'process_cluster' repeatedly. The transferred identifications
are stored in the 'identification_references' member variable.

@author Johannes Griss
"""

from . import common
from .. import objects


class IdTransferer(common.AbstractAnalyser):
    def __init__(self, add_to_identified=False, add_to_unidentified=True, include_all_identified=False):
        """
        Creates a default IdTransferer object.
        :param add_to_identified: If set identifications are added to identified spectra.
        :param add_to_unidentified: If set identifications are added to unidentified spectra.
        :param include_all_identified: If set identified spectra that are not part of reliable clusters
         are returned as well. Additionally, if add_to_identified is set to false and
         include_all_identified is set to true, the original identifications are returned unchanged.
        :return:
        """
        super().__init__()
        self.identification_references = list()
        self.add_to_identified = add_to_identified
        self.add_to_unidentified = add_to_unidentified
        self.include_all_identified = include_all_identified

    def process_cluster(self, cluster):
        """
        Transfers ids to spectra based on the cluster's properties
        :param cluster: The cluster to process
        :return:
        """
        # this only works on identified clusters
        if cluster.identified_spectra < 1:
            return

        # test whether the cluster should be processed at all
        if self._ignore_cluster(cluster):
            if self.include_all_identified:
                # keep all identified spectra
                for spectrum in cluster.get_spectra():
                    if spectrum.is_identified():
                        # just use the existing ids
                        self.identification_references.append(
                            IdentificationReference(
                                spectrum.get_filename(),
                                spectrum.get_id(),
                                spectrum.get_clean_sequence_psms()))
            # do not process this cluster any further
            return

        # get the primary identification object
        main_psms = IdTransferer.extract_main_cluster_psms(cluster)

        # add it to all spectra
        for spectrum in cluster.get_spectra():
            # save the original identification if keep all identified is set
            if spectrum.is_identified() and self.include_all_identified and not self.add_to_identified:
                self.identification_references.append(
                    IdentificationReference(
                        spectrum.get_filename(),
                        spectrum.get_id(),
                        spectrum.get_clean_sequence_psms()))
                continue

            # make sure the identification should be added to the spectrum
            if not self.add_to_identified and spectrum.is_identified():
                continue
            if not self.add_to_unidentified and not spectrum.is_identified():
                continue

            # create and add the new identification reference
            self.identification_references.append(
                IdentificationReference(spectrum.get_filename(), spectrum.get_id(), main_psms))

    @staticmethod
    def extract_main_cluster_psms(cluster):
        """
        Extracts the most common identification(s) from a cluster
        and returns it as a PSM object.
        :param cluster: Cluster to extract the most common PSM from.
        :return: Most common PSMs as a list
        """
        main_psms = list()

        for max_sequence in cluster.max_sequences:
            # get the modifications for this sequence
            ptms = list()

            for spectrum in cluster.get_spectra():
                for psm in spectrum.psms:
                    if psm.sequence == max_sequence and len(psm.ptms) > len(ptms):
                        ptms = psm.ptms

            main_psms.append(objects.PSM(max_sequence, ptms))

        return main_psms


class IdentificationReference:
    """Contains the main identification data."""

    def __init__(self, filename, spec_id, psms):
        """
        Creates a new instance of the identification reference.
        :param filename: Original peak list filename.
        :param spec_id: The spectrum's id within this file.
        :param psms: A list of PSM objects
        :return:
        """
        self.filename = filename
        self.spec_id = spec_id
        self.psms = psms