"""
This analyser extracts the number of spectra per
sample and cluster to output a table containing
the samples as columns and the clusters as rows.
"""

from . import common
import pandas


class ClusterAsFeatures(common.AbstractAnalyser):
    """
    Extracts the number of spectra per
    sample and cluster to output a table containing
    the samples as columns and the clusters as rows.

    :ivar features: A list of dicts. Each dict represents
                    one cluster. The dict contains the sample
                    name as key and the number of spectra as value.

    TODO: add description about samples
    """
    def __init__(self, sample_name_extractor=None):
        """
        Initialised a new ClusterAsFeatures analyser.
        :param sample_name_extractor: A function that takes
        the spectrum's title as parameter and returns the ie.
        sample name. If set to None the default function is used
        where everything before the first "." is being returned.
        :return:
        """
        super().__init__()

        if sample_name_extractor is None:
            sample_name_extractor = ClusterAsFeatures.extractBasicSampleName

        self.sample_name_extractor = sample_name_extractor
        self.features = list()
        self.cluster_ids = list()
        self.samples = set()

    @staticmethod
    def extractBasicSampleName(spec_ref):
        """
        Extracts the sample name by returning everything before
        the first "." from the title (often used by ProteoWized
        converted files)

        :param spec_ref: The spectrum object.
        :return: The sample name
        """
        spectrum_title = spec_ref.get_title()
        index = spectrum_title.find(".")

        if index < 0:
            return spectrum_title

        return spectrum_title[index + 1:]

    def process_cluster(self, cluster):
        """
        Extracts how many spectrum per sample were
        observed.

        :param cluster: The cluster to process
        :return:
        """
        if self._ignore_cluster(cluster):
            return

        spec_per_sample = dict()

        for spec_ref in cluster.get_spectra():
            sample_id = self.sample_name_extractor(spec_ref)

            # collect all unique samples
            self.samples.add(sample_id)

            if sample_id in spec_per_sample:
                spec_per_sample[sample_id] += 1
            else:
                spec_per_sample[sample_id] = 1

        self.features.append(spec_per_sample)
        self.cluster_ids.append(cluster.id)

    def get_result(self):
        """
        Return the result as a pandas DataFrame.

        :return: A numpy array representing the merged result.
        """
        sample_list = list(self.samples)
        result = pandas.DataFrame(columns=sample_list, index=self.cluster_ids)

        for i in range(0, len(self.cluster_ids)):
            cluster_id = self.cluster_ids[i]
            spectra_per_sample = self.features[i]
            for sample_id in sample_list:
                result.loc[cluster_id, sample_id] = spectra_per_sample.get(sample_id, 0)

        return result



