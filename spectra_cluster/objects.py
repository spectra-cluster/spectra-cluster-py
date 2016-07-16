# -----------------------------
# This contains the data types
# used by the clustering parser.
# ------------------------------

import re


class Cluster:
    """
    Represents a cluster in a .clustering output file.
    """

    # TODO: add helper functions: (i/l) ratio, ident spec count, unid spec count, has ident spec

    def __init__(self, id, precursor_mz, consensus_mz, consensus_intens, spectra):
        """
        Creates a new cluster object
        :param id: The cluster's id
        :param precursor_mz: The cluster's average precursor m/z
        :param consensus_mz: A list of doubles holding the consensus spectrum's m/z values
        :param consensus_intens: A list of doubles holding the consensus spectrum's intensity values
        :param spectra: A list of spectra associated with the cluster
        :return:
        """
        self.id = id
        self.precursor_mz = precursor_mz
        self.consensus_mz = consensus_mz
        self.consensus_intes = consensus_intens
        self.spectra = spectra


class PSM:
    """
    A peptide-spectrum-match reference.
    """

    def __init__(self, title, sequences, precursor_mz, charge, taxids, ptms):
        """
        Creates a new PSM reference.
        :param title: The PSM's title.
        :param sequences: Sequences associated with this spectrum.
        :param precursor_mz: Measured precursor m/z
        :param charge: Charge state
        :param taxids: Taxids of the experiments in which the spectrum was observed
        :param ptms: A list of PTMs
        :return:
        """
        self.title = title
        self.sequences = sequences
        self.precursor_mz = precursor_mz
        self.charge = charge
        self.taxids = taxids
        self.ptms = ptms

    def get_filename(self):
        """
        The originally filename can optionally be encoded
        in the title string. If present this filename is returned
        otherwise None
        :return: Original filename or None if not present
        """
        if "#file=" not in self.title:
            return None

        start = self.title.find("#file=")
        end = self.title.find("#id=")

        if end == -1:
            end = len(self.title)

        return self.title[start + 6:end]

    def get_id(self):
        """
        The spectrum's id can optionally be encoded
        in the title string. If present this id is
        returned, otherwise None.
        :return: Original spectrum id or None if not present
        """
        if "#id=" not in self.title:
            return None

        start = self.title.find("#id=")
        end = self.title.find("#title=")

        if end == -1:
            end = len(self.title)

        return self.title[start + 4:end]

    def get_title(self):
        """
        Optionally the spectrum's filename and id can be encoded
        in the title string. If this is the cases, the original
        title is extracted from the string. If no fields were
        encoded, the whole title string is returned. Therefore, this
        function should always be used if the reader expects
        to access the original spectrum's title.
        :return: The original spectrum's title
        """
        if "#title=" not in self.title:
            return self.title

        start = self.title.find("#title=")
        return self.title[start + 7:]

    def get_clean_sequences(self):
        """
        Returns the identified sequences without any additional
        characters and only using high-caps.
        :return: Identified sequences
        """
        clean_sequences = list()

        for sequence in self.sequences:
            clean_sequences.append(
                re.sub(r"[^A-Z]", "", sequence.upper()))

        return clean_sequences


class PTM:
    """
    Defines a post-translational modification within a peptide
    """

    def __init__(self, position, accession):
        """
        Creates a new PTM object
        :param position: 1-based position within the peptide (0 for terminus)
        :param accession: MOD accession of the modification.
        :return:
        """
        self.position = position
        self.accession = accession
