"""
This exporter writes the consensus spectra of the
clusters into an MGF formatted file.
"""

from .. import common


class MgfExporter(common.AbstractAnalyser):
    """
    Converts the clusters' consensus spectra into MGF format.
    """
    def __init__(self, result_file):
        """
        Initialises a new MgfExporter object.

        :param result_file: File object to write to.
        """
        super().__init__()

        self.result_file = result_file

    def process_cluster(self, cluster):
        """
        Convert the cluster's consensus spectrum
        into MGF format.

        :param cluster: The cluster to process
        """
        if self._ignore_cluster(cluster):
            return

        # make sure the peak list is valid
        if len(cluster.consensus_mz) != len(cluster.consensus_intens):
            raise Exception("Cluster " + cluster.id + " contains different number of m/z and"
                                                      "intensity values for the consensus spectrum")

        # build the title string
        title_string = "TITLE=" + cluster.id
        if cluster.identified_spectra > 0:
            title_string += ",sequence=" + ",".join(cluster.max_sequences)
        else:
            title_string += ",sequence=UNIDENTIFIED"

        # write the spectrum
        lines = list(["BEGIN IONS", title_string])

        lines.append("PEPMASS=" + str(round(cluster.precursor_mz, 4)))
        lines.append("CHARGE=" + str(cluster.charge))

        if cluster.identified_spectra > 0:
            lines.append("SEQUENCE=" + ",".join(cluster.max_sequences))

        # add the peak list
        for i in range(0, len(cluster.consensus_mz)):
            lines.append(str(round(cluster.consensus_mz[i], ndigits=4)) + " " +
                         str(round(cluster.consensus_intens[i], ndigits=4)))

        lines.append("END IONS\n\n")

        # write the spectrum
        self.result_file.write("\n".join(lines))
