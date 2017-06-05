from . import common
from .. import objects


class IdTransferer(common.AbstractAnalyser):
    """
    The IdTransferer analyser transfers identifications to spectra part of a cluster.

    The analysis is run by calling 'process_cluster' repeatedly. The transferred identifications
    are stored in the 'identification_references' member variable.

    :ivar identification_references: A list of IdentificationReferences
    """
    def __init__(self, add_to_identified=False, add_to_unidentified=True, include_all_identified=False):
        """
        Creates a default IdTransferer object.

        :param add_to_identified: If set identifications are added to identified spectra.
        :param add_to_unidentified: If set identifications are added to unidentified spectra.
        :param include_all_identified: If set identified spectra that are not part of reliable clusters
         are returned as well. Additionally, if add_to_identified is set to false and
         include_all_identified is set to true, the original identifications are returned unchanged.

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
                                spectrum.get_clean_sequence_psms(),
                                False,
                                spectrum=spectrum))
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
                        spectrum.get_clean_sequence_psms(),
                        False,
                        spectrum=spectrum))
                continue

            # make sure the identification should be added to the spectrum
            if not self.add_to_identified and spectrum.is_identified():
                continue
            if not self.add_to_unidentified and not spectrum.is_identified():
                continue

            # check whether the identification details are changed
            changed_identification = False
            if not spectrum.is_identified():
                changed_identification = True
            else:
                sequences = [s for s in spectrum.get_clean_sequences()]
                if sequences[0] != main_psms[0].sequence:
                    changed_identification = True

            # create and add the new identification reference
            self.identification_references.append(
                IdentificationReference(spectrum.get_filename(),
                                        spectrum.get_id(), main_psms,
                                        changed_identification,
                                        spectrum=spectrum))

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
    """
    Contains the main identification data.

    :ivar filename: The original peak list filename
    :ivar spec_id: The spectrum's id within the source file
    :ivar psms: A list of PSM objects
    :ivar changed_through_clustering: Logical indicating whether
          the identification details were changed through the
          clustering.
    :ivar spectrum: The complete spectrum object
    """

    def __init__(self, filename, spec_id, psms, changed_through_clustering, spectrum):
        """
        Creates a new instance of the identification reference.

        :param filename: Original peak list filename.
        :param spec_id: The spectrum's id within this file.
        :param psms: A list of PSM objects
        :param changed_through_clustering: Logical indicating whether
          the identification details were changed through the
        """
        self.filename = filename
        self.spec_id = spec_id
        self.psms = psms
        self.changed_through_clustering = changed_through_clustering
        self.spectrum = spectrum