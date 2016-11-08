"""
This analyser extracts the number of spectra per
sample and cluster to output a table containing
the samples as columns and the clusters as rows.
"""

from . import common
import operator  


class ClusterListsComparer(common.AbstractAnalyser):
    """
    This tool compare two cluster lists 
    and give the statistics between them.

    Result
    ------
    The results are stored in ::tableList:: as a list of
    tables. Each table represents a statistics information. 

    TODO: 
    """
    def __init__(self):
        """
        Initialised a new ClusterListsComparer analyser.

        :return:
        """
        super().__init__()
        self.cluster_lists = [[],[]] 
        self.tables = dict()
        self.samples = set()
        self.dealing_list = 0
        self.sorted_spectra_dict = dict()

        self.shared_spec_num = 0
        self.similarity_dis = dict()
        self.cluster_added = 0
    def process_cluster(self, cluster):
        if self._ignore_cluster(cluster):
            return
        self.cluster_added += 1
        self.cluster_lists[self.dealing_list].append(cluster)
        spectra = list(cluster.get_spectra())

        def mixed_order( spec ):
            return ( spec.get_filename(), spec.get_id())
#            return ( spec.get_filename(), spec.get_id(), spec.get_clean_sequences())
#        print("\n\n--------")
#        print("befor-sorted")
#        for spec in spectra:
#            print(spec.get_filename() + spec.get_id())
        spectra.sort(key=mixed_order) 
#        print("after-sorted")
#        for spec in spectra:
#            print(spec.get_filename() + spec.get_id())
        self.sorted_spectra_dict[cluster.id] = spectra

        
    def calculate_similarity(self, cluster0, cluster1):

        ###########################################
        def compare_spectrum(spectrum0, spectrum1):
            filename0 = spectrum0.get_filename() 
            filename1 = spectrum1.get_filename() 
            id0 = spectrum0.get_id()
            id1 = spectrum1.get_id()

            if(filename0 < filename1): 
                return -1
            elif(filename0 > filename1):   
                return 1
            elif(id0 < id1):
                return -1
            elif(id0 > id1): 
                return 1 
            else:
                return 0
        ###########################################

        spectra0 = self.sorted_spectra_dict[cluster0.id]
        spectra1 = self.sorted_spectra_dict[cluster1.id]

        (n,i,j) = (0,0,0)
        while(i<len(spectra0) and j<len(spectra1)):
            comp_score = compare_spectrum(spectra0[i], spectra1[j])
            if(comp_score < 0):
                i += 1
            elif(comp_score > 0):
                j += 1
            else:    
                n += 1
                i += 1
                j += 1
        similarity_score = 0.5 * (n/len(spectra0) + n/len(spectra1))
        self.shared_spec_num += n
        return similarity_score


    def compare(self):
        (i,j) = (0,0)
        for cluster0 in self.cluster_lists[0]:
            i += 1
            j = 0
            for cluster1 in self.cluster_lists[1]:
                j += 1
                similarity = self.calculate_similarity(cluster0, cluster1)
                self.similarity_dis[str(int(similarity*10))] = self.similarity_dis.get(str(int(similarity*10)),0) + 1
#                if(similarity == 1):
#                    print( str(i) +":"+ str(j) + "-"+cluster0.id +"-" +cluster1.id + "--"+str(self.similarity_dis[str(int(similarity*10))]))

    def output_debug_info(self):
        print(len(self.cluster_lists[0]))
        print(len(self.cluster_lists[1]))
        print(self.similarity_dis)
        print(self.shared_spec_num)
        print(self.cluster_added )

        
