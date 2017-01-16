"""
This analyser extracts the number of spectra per
sample and cluster to output a table containing
the samples as columns and the clusters as rows.
"""

from . import common
import tempfile
import shutil
import os


class ClusterAsFeatures(common.AbstractAnalyser):
    """
    Extracts the number of spectra per
    sample and cluster and writes the result directly
    to a file object.

    Since the number of samples within the clustering
    result is not known at the beginning, you have to
    use the function "add_resultfile_header" to add a
    header to the result file.
    """
    def __init__(self, result_file, sample_name_extractor=None):
        """
        Initialised a new ClusterAsFeatures analyser.

        :param result_file: A file object that will be used to
                            write the resulting table to. This is
                            necessary since this result data will
                            generally be too large to keep in memory.
        :param sample_name_extractor: A function that takes the spectrum's
                                      title as parameter and returns the ie.
                                      ample name. If set to None the default
                                      function is used where everything before
                                      the first "." is being returned.
        """
        super().__init__()

        if sample_name_extractor is None:
            sample_name_extractor = ClusterAsFeatures.extract_basic_sample_name

        self.sample_name_extractor = sample_name_extractor
        self.result_file = result_file
        self.sample_ids = list()

    @staticmethod
    def extract_basic_sample_name(spec_ref):
        """
        Extracts the sample name by returning everything before
        the first "." from the title (often used by ProteoWized
        converted files) or, if available, by returning the original
        filename (without path information).

        :param spec_ref: The spectrum object.
        :return: The sample name
        """
        # use the filename if available
        spectrum_filename = spec_ref.get_filename()

        if spectrum_filename is not None:
            return os.path.basename(spectrum_filename)

        # otherwise use the spectrum title
        spectrum_title = spec_ref.get_title()
        index = spectrum_title.find(".")

        if index < 0:
            return spectrum_title

        return spectrum_title[:index]

    def process_cluster(self, cluster):
        """
        Extracts how many spectrum per sample were
        observed.

        :param cluster: The cluster to process
        :return:
        """
        if self._ignore_cluster(cluster):
            return

        # count the number of spectra per sample
        spec_per_sample = dict()

        for spec_ref in cluster.get_spectra():
            sample_id = self.sample_name_extractor(spec_ref)

            if sample_id in spec_per_sample:
                spec_per_sample[sample_id] += 1
            else:
                spec_per_sample[sample_id] = 1

        # add all samples that haven't been added yet
        for sample_id in spec_per_sample.keys():
            if sample_id not in self.sample_ids:
                self.sample_ids.append(sample_id)

        # write the table - first column is always the cluster id
        fields = [cluster.id]
        for sample_id in self.sample_ids:
            fields.append(str(spec_per_sample.get(sample_id, 0)))

        result_line = "\t".join(fields)
        self.result_file.write(result_line + "\n")

    def add_resultfile_header(self, file_path):
        """
        Adds the header line to the result file that
        was used to write the results during the analysis
        to.

        :param file_path: Path ot the file where the results
                          are stored.
        """
        tmp = tempfile.TemporaryFile(mode="w+")

        # write the header to the temporary file
        fields = ["cluster_id"]
        fields += self.sample_ids

        field_string = "\t".join(fields) + "\n"

        tmp.write(field_string)

        # copy the result
        with open(file_path, "r") as IN:
            for line in IN:
                tmp.write(line)

        tmp.seek(0)

        with open(file_path, "w") as OUT:
            shutil.copyfileobj(tmp, OUT)

#    def get_result(self):
#        """
#        Return the result as a pandas DataFrame.
#
#        :return: A numpy array representing the merged result.
#        """
#        sample_list = list(self.samples)
#        result = DataFrame(columns=sample_list, index=self.cluster_ids)
#
#        for i in range(0, len(self.cluster_ids)):
#            cluster_id = self.cluster_ids[i]
#            spectra_per_sample = self.features[i]
#            for sample_id in sample_list:
#                result.loc[cluster_id, sample_id] = spectra_per_sample.get(sample_id, 0)
#
#        return result



