# -----------------------------
# This contains the data types
# used by the clustering parser.
# ------------------------------

import re


class Cluster:
    """
    Represents a cluster in a .clustering output file.
    """
    def __init__(self, cluster_id, precursor_mz, consensus_mz, consensus_intens, spectra):
        """ Creates a new cluster object

        :param cluster_id: The cluster's id
        :param precursor_mz: The cluster's average precursor m/z
        :param consensus_mz: A list of doubles holding the consensus spectrum's m/z values
        :param consensus_intens: A list of doubles holding the consensus spectrum's intensity values
        :param spectra: A set of spectra associated with the cluster
        :return:
        """
        self.id = cluster_id
        self.precursor_mz = precursor_mz
        self.consensus_mz = consensus_mz
        self.consensus_intens = consensus_intens
        self._spectra = set(spectra)

        self._update_properties()

    def _update_properties(self):
        """
        This function calculates additional properties such as most common
        sequence, max I/L ignorant ratio, number of identified spectra etc.
        It should always be called if the spectra object is changed.

        :return:
        """
        self.n_spectra = len(self._spectra)
        self.identified_spectra = 0
        self.unidentified_spectra = 0

        for spectrum in self._spectra:
            if spectrum.is_identified():
                self.identified_spectra += 1
            else:
                self.unidentified_spectra += 1

        if self.identified_spectra > 0:
            # calculate ratios
            sequence_counts = Cluster.calculate_sequence_counts(self._spectra, False)

            self.sequence_ratios = dict()
            for sequence in sequence_counts.keys():
                self.sequence_ratios[sequence] = sequence_counts[sequence] / self.identified_spectra

            self.max_ratio = max(self.sequence_ratios.values())

            self.max_sequences = tuple([sequence for sequence in self.sequence_ratios.keys()
                                        if self.sequence_ratios[sequence] == self.max_ratio])

            # calculate max I/L ratio
            sequence_counts_il = Cluster.calculate_sequence_counts(self._spectra, True)

            sequence_ratios_il = dict()
            for sequence in sequence_counts_il.keys():
                sequence_ratios_il[sequence] = sequence_counts_il[sequence] / self.identified_spectra

            self.max_il_ratio = max(sequence_ratios_il.values())
        else:
            # set to default values for unidentified clusters
            self.sequence_ratios = dict()
            self.max_ratio = None
            self.max_sequences = tuple()
            self.max_il_ratio = None

    def get_spectra(self):
        """
        Returns the stored spectra in a tuple. These object should
        not be changed. Otherwise, the cluster's statistics may no
        longer be accurate.

        :return: A tuple containing the cluster's spectra
        """
        return tuple(self._spectra)

    def set_spectra(self, new_spectra):
        """
        Updates the cluster's stored spectra

        :param new_spectra: A list of PSM objects.
        """
        self._spectra = list(new_spectra)
        self._update_properties()

    @staticmethod
    def calculate_sequence_counts(spectra, ignore_i_l=False):
        """
        Calculates the sequence counts based on
        the passed spectra. PTMs are ignored for this
        assessment

        :param spectra: The spectra to derive the sequence counts from.
        :param ignore_i_l: If set I and L are treated as equivalent. If set all I are
          replace by L and the sequences in the returned map may not correspond to the
          originally identified sequences.
        :return: A dict with a sequence as key and the number of occurrences as value.
        """
        sequence_counts = dict()

        for spectrum in spectra:
            if not spectrum.is_identified():
                continue

            for sequence in spectrum.get_clean_sequences():
                if ignore_i_l:
                    processed_sequence = sequence.replace("I", "L")
                else:
                    processed_sequence = sequence

                if processed_sequence not in sequence_counts:
                    sequence_counts[processed_sequence] = 1
                else:
                    sequence_counts[processed_sequence] += 1

        return sequence_counts


class Spectrum:
    """
    A spectrum reference.
    """

    def __init__(self, title, precursor_mz, charge, taxids, psms):
        """
        Creates a new Spectrum reference.

        :param title: The spectrum's title.
        :param precursor_mz: Measured precursor m/z
        :param charge: Charge state
        :param taxids: Set of taxids of the experiments in which the spectrum was observed
        :param psms: A set of psms associated with the spectrum. If None is passed
        an empty set is created.
        :return:
        """
        self.title = title
        self.precursor_mz = precursor_mz
        self.charge = charge
        self.taxids = frozenset(taxids)

        if psms is None:
            self.psms = frozenset()
        else:
            self.psms = frozenset(psms)

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

    def get_clean_sequence_psms(self):
        """
        Returns all PSMs with all special characters
        removed from the sequences.

        :return: A tuple of PSMs
        """
        clean_psms = list()

        for psm in self.psms:
            clean_sequence = re.sub(r"[^A-Z]", "", psm.sequence.upper())
            clean_psms.append(PSM(sequence=clean_sequence, ptms=psm.ptms))

        return tuple(clean_psms)

    def get_clean_sequences(self):
        """
        Returns the identified sequences without any additional
        characters and only using high-caps.

        :return: Identified sequences
        """
        clean_sequences = set()

        for psm in self.psms:
            clean_sequences.add(
                re.sub(r"[^A-Z]", "", psm.sequence.upper()))

        return clean_sequences

    def is_identified(self):
        """
        Checks whether the spectrum was identified.

        :return: boolean
        """
        if self.psms is None:
            return False

        return len(self.psms) > 0

    def __eq__(self, other):
        """
        Tests whether this Spectrum class is describing the same spectrum as the
        passed other instance.

        :param other: Object to compare to
        :return: Boolean indicating whether both object are equal.
        """
        if type(self) != type(other):
            return False

        if self.title != other.title:
            return False
        if self.precursor_mz != other.precursor_mz:
            return False
        if self.charge != other.charge:
            return False
        if self.taxids != other.taxids:
            return False
        if self.psms != other.psms:
            return False

        return True

    def __str__(self):
        return "<Spectrum @ " + str(self.precursor_mz) + " m/z " + str(self.charge) + \
               "+ with " + str(len(self.psms)) + " PSMs>"

    def __repr__(self):
        return repr(self.title) + \
               repr(self.precursor_mz) + \
               repr(self.charge) + \
               repr(self.taxids) + \
               repr(self.psms)

    def __hash__(self):
        return hash(repr(self))


class PSM:
    """
    Defines a peptide-spectrum-match
    """
    def __init__(self, sequence, ptms):
        """
        Creates a new PSM object.

        :param sequence: The sequence associated with the PSM.
        :param ptms: A set of PTMs
        :return:
        """
        self.sequence = sequence
        self.ptms = frozenset(ptms)

    def __eq__(self, other):
        """
        PSMs are defined based on the sequence and the
        PTMs.

        :param other: The class to compare to
        :return: Boolean indicating whether both objects represent the same PSM.
        """
        if type(self) != type(other):
            return False

        if self.sequence != other.sequence:
            return False

        if self.ptms != other.ptms:
            return False

        return True

    def __str__(self):
        sequence_list = list()
        for index in range(0, len(self.sequence)):
            for ptm in self.ptms:
                if ptm.position == index:
                    sequence_list.append("[" + ptm.accession + "]")

            sequence_list.append(self.sequence[index])

        # check the last termini
        for ptm in self.ptms:
            if ptm.position == len(self.sequence):
                sequence_list.append("[" + ptm.accession + "]")

        return "".join(sequence_list)

    def __repr__(self):
        return repr(self.sequence) + repr(self.ptms)

    def __hash__(self):
        return hash(repr(self))


class PTM:
    """
    Defines a post-translational modification within a peptide

    :ivar position: The PTM's position within the peptide string
    :ivar accession: The PTM's accession in UNIMOD (if starting with
                     "MOD:"). This may also represent a PSI entry in
                     the format [PSI-MS, MS:1001524, fragment neutral loss, 63.998283]
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

    def __eq__(self, other):
        """
        PTMs are defined based on the position and the
        accession.

        :param other: The class to compare to
        :return: Boolean indicating whether both objects represent the same PTM.
        """
        if type(self) != type(other):
            return False

        if self.accession != other.accession:
            return False

        if self.position != other.position:
            return False

        return True

    def __repr__(self):
        return repr(self.position) + repr(self.accession)

    def __hash__(self):
        return hash(repr(self))

    def __str__(self):
        return str(self.position) + "-" + self.accession
