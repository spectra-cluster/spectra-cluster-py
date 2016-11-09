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
        self.min_size = 2 #set default 2
        self.cluster_lists = [[],[]] 
        self.cluster_spectra_num = [0,0] 
        self.distr_cluster_size = [{},{}]
        self.tables = list()
        self.dealing_list = 0
        self.sorted_spectra_dict = dict()
        self.stars = dict()
        self.starlets = dict()

        self.shared_spec_num = 0
        self.similarity_dis = dict()
        self.stars_length = 0
        self.starlets_length = 0
        self.star_divide_factor_dist = dict()
        self.standalone_star_num = 0


    def process_cluster(self, cluster):
        if self._ignore_cluster(cluster):
            return
        self.cluster_lists[self.dealing_list].append(cluster)
        spectra = list(cluster.get_spectra())

        def mixed_order( spec ):
            return ( spec.get_filename(), spec.get_id())
        spectra.sort(key=mixed_order) 
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
        return (n,similarity_score)


    def compare(self):
        similarity_dis = dict()
        i = 0
        len0 = len(self.cluster_lists[0])
        len1 = len(self.cluster_lists[1])
        self.stars_length = len0 if len0<=len1 else len1
        self.starlets_length = len1 if len0>len1 else len0

        for cluster in self.cluster_lists[0]:
            len_spectra =  len(cluster.get_spectra())
            if len0 <= len1:
                star = ClusterNode(cluster.id, len_spectra) 
                self.stars[cluster.id] = star
            else:
                starlet = ClusterNode(cluster.id, len_spectra) 
                self.starlets[cluster.id] = starlet

            self.cluster_spectra_num[0] += len_spectra
            self.distr_cluster_size[0][len_spectra] = self.distr_cluster_size[0].get(len_spectra,0) + 1

        for cluster in self.cluster_lists[1]:
            len_spectra =  len(cluster.get_spectra())
            if len0 <= len1:
                starlet = ClusterNode(cluster.id, len_spectra) 
                self.starlets[cluster.id] = starlet
            else:
                star = ClusterNode(cluster.id, len_spectra) 
                self.stars[cluster.id] = star

            self.cluster_spectra_num[1] += len_spectra
            self.distr_cluster_size[1][len_spectra] = self.distr_cluster_size[1].get(len_spectra,0) + 1

        for cluster0 in self.cluster_lists[0]:
            i += 1
            j = 0
            for cluster1 in self.cluster_lists[1]:
                j += 1
                (shared_spec_num, similarity) = self.calculate_similarity(cluster0, cluster1)
                if similarity == 0:
                    continue
                self.similarity_dis[str(int(similarity*10))] = self.similarity_dis.get(str(int(similarity*10)),0) + 1
                self.shared_spec_num += shared_spec_num
                if len0 <= len1:
                    self.stars[cluster0.id].add_nb_node(cluster1.id, similarity, shared_spec_num)
                    self.starlets[cluster1.id].add_nb_node(cluster0.id, similarity, shared_spec_num)
                else:
                    self.stars[cluster1.id].add_nb_node(cluster0.id, similarity, shared_spec_num)
                    self.starlets[cluster0.id].add_nb_node(cluster1.id, similarity,shared_spec_num)

        self.shared_spec_num /= 2 

    
    def caculate_network_statistics(self):
        self.standalone_star_num += 1 
        for key in self.stars.keys():
            star = self.stars[key]
            if star.nb_num == 0 :
                self.standalone_star_num += 1 
                continue
            divide_factor = star.nb_num + star.spec_num - star.shared_spec_num
            self.star_divide_factor_dist[divide_factor] = self.star_divide_factor_dist.get(divide_factor,0) + 1
            if star.spec_num < star.shared_spec_num:
                print("!!!!!!!!!!!!Wrong here, total spectra No is less than Shared Spectra with starlets")
                print("with star " + star.id + "  " + str(star.spec_num) + "is less than" + str(star.shared_spec_num))
            if star.spec_num > star.shared_spec_num:
                star_lost_spec_numt += star.spec_num - star.shared_spec_num
        



    def prepare_statistics(self):
        tabs = "\t\t" 
        #statistics of clustering files
        percentage0 = 100.0 * float(self.shared_spec_num)/float(self.cluster_spectra_num[0])
        percentage1 = 100.0 * float(self.shared_spec_num)/float(self.cluster_spectra_num[1])
        rows = "name %s number %s description\n" %(tabs, tabs)
        rows += "cluster No %s%s %d %s in file0\n" % (tabs, tabs, len(self.cluster_lists[0]), tabs)
        rows += "cluster No %s%s %d %s in file1\n" % (tabs, tabs, len(self.cluster_lists[1]), tabs)
        rows += "identical cluster No %s%s %s %s between them\n" % (tabs, tabs, self.similarity_dis['10'], tabs)
        rows += "spectrum No %s %d %s in file0\n" % (tabs, self.cluster_spectra_num[0], tabs)
        rows += "spectrum No %s %d %s in file1\n" % (tabs,self.cluster_spectra_num[1], tabs)
        rows += "shared spectrum No %s %d %s between them\n" % (tabs, self.shared_spec_num, tabs)
        rows += "shared spectrum percent %s %.2f %s in file0\n" % (tabs, percentage0, tabs)
        rows += "shared spectrum percent %s %.2f %s in file1\n" % (tabs, percentage1, tabs)
        self.tables.append(('statistics of files', rows))

        #distribution of cluster size in file0
        rows = "cluster size %s No. %s percentage %s accumulate pecent\n" % (tabs, tabs, tabs)
        len0 = len(self.cluster_lists[0])
        accumulate_num = 0
        for key in self.distr_cluster_size[0].keys():
            value = self.distr_cluster_size[0][key]
            accumulate_num += value
            percent = 100 * value/len0
            accum_percent = 100 * accumulate_num/len0
            rows += "%d %s %d %s %.2f %s %.2f\n" % (key, tabs, value, tabs, percent, tabs,accum_percent)
        self.tables.append(('distribution of cluster size in file0', rows))
        
        rows = "cluster size %s No. %s percentage %s accumulate pecent\n" % (tabs, tabs, tabs)
        len1 = len(self.cluster_lists[1])
        accumulate_num = 0
        for key in self.distr_cluster_size[1].keys():
            value = self.distr_cluster_size[1][key]
            accumulate_num += value
            percent = 100 * value/len1
            accum_percent = 100 * accumulate_num/len0
            rows += "%d %s %d %s %.2f %s %.2f\n" % (key, tabs, value, tabs, percent, tabs,accum_percent)
        self.tables.append(('distribution of cluster size in file1', rows) )

        #distribution of similarity
        rows = "similarity score %s pairs of clusters %s percentage(stars) %s percentage(starlets)\n" % (tabs, tabs, tabs)
        for key in self.similarity_dis.keys():
            value = self.similarity_dis[key]
            rows += "%d %s %d %s %.2f %s %.2f\n" % (int(key), tabs, int(value), tabs, 100.0*value/self.stars_length, tabs, 100.0*value/self.starlets_length)
        self.tables.append(('distribution of similarity (identical = 10)', rows))

        #distribution of star divide factors
        rows = "divide factor %s No\n" % (tabs) 
        for key in self.star_divide_factor_dist.keys():
            value = self.star_divide_factor_dist[key]
            rows += "%d %s %d \n" %(int(key), tabs, int(value) )
        self.tables.append(('distribution of star divide factors', rows))

    def output_debug_info(self):
        print(len(self.cluster_lists[0]))
        print(len(self.cluster_lists[1]))
#        print(self.tables['similarity_dis'])
        print(self.shared_spec_num)

#        print(self.stars)
#        print(self.starlets)

class ClusterNode:  
    """
    Simple class to represent the Clusters as nodes in the connected network.
    """

    def __init__(self, cluster_id, spec_num):
        self.id = cluster_id
        self.spec_num = spec_num
        self.shared_spec_num = 0 

        self.nb_num = 0
        self.connectedNodes = list()

    def add_nb_node(self,nb_id, similarity, shared_spec_num):
        self.connectedNodes.append((nb_id,similarity, shared_spec_num))
        self.nb_num += 1
        self.shared_spec_num += shared_spec_num 
