"""
This analyser compare the diffrence between two
cluster lists, give the statsistics tables 
about the  differences.
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
        self.min_size = 2 # set default minium size 2
        self.file_index = 0
       
        # intermediate data structures
        self.cluster_lists = [[],[]] 
        self.cluster_spectra_num = [0,0] 
        self.sorted_spectra_dict = dict()
        self.stars = dict()
        self.starlets = dict()

        # statisitic resluts
        self.shared_spec_num = 0
        self.stars_length = 0
        self.starlets_length = 0
        self.standalone_star_num = 0
        self.standalone_starlet_num = 0
        self.similarity_dist = dict()
        self.star_divide_factor_dist = dict()
        self.starlet_divide_factor_dist = dict()
        self.cluster_size_dist = [{},{}]
        
        # output tables
        self.tables = list()

    def process_cluster(self, cluster):
        """
        Add the clusters into two cluster lists, based on 
        the files they come from. 
        
        The spectra list in each cluster has been sorted
        based on the "filename" and "index(spec_id)", to 
        reduce the comparing cost. Cluster with sorted 
        spectra is stored in self.sorted_spectra_dict.
        
        :param cluster: the cluster to process
        :return:
        """
        if self._ignore_cluster(cluster):
            return

        self.cluster_lists[self.file_index].append(cluster)
        spectra = list(cluster.get_spectra())

        def mixed_order( spec ): # for sorting
            return ( spec.get_filename(), spec.get_id())

        spectra.sort(key=mixed_order) 
        self.sorted_spectra_dict[cluster.id] = spectra

        
    def calculate_similarity(self, cluster0, cluster1):
        """
        Calculate the similarity between to cluster.
        
        :param cluster0 cluster1: the pair of clusters
        to calculate

        :return(n,similarity_score): n is the shared
        spectra number tetween them, similarity_score
        is explained by name. 
        """

        def compare_spectrum(spectrum0, spectrum1): 
            """
            Compare a pair of spectra to decide the
            order. 
            :param: pair of spectra
            :return: 0 equal, -1 spectrum0 is less,
            1, spectrum0 is bigger.
            """
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
        # end of compare_spectrum function

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
        """
        Compare two cluster lists, get the distribution of the similarities.
        And build the connection network between clusters.
        """
        similarity_dist = dict()
        len0 = len(self.cluster_lists[0])
        len1 = len(self.cluster_lists[1])
        longer_index = 0 if len0 >= len1 else 1
        shorter_index = 1 if len1 > len0 else 0
        self.stars_length = len0 if len0<=len1 else len1
        self.starlets_length = len1 if len0>len1 else len0

        # build the noeds for cluster list0, and get the
        # distribution of cluster size.
        for cluster in self.cluster_lists[0]:
            len_spectra =  len(cluster.get_spectra())
            if len0 <= len1:
                star = ClusterNode(cluster.id, len_spectra) 
                self.stars[cluster.id] = star
            else:
                starlet = ClusterNode(cluster.id, len_spectra) 
                self.starlets[cluster.id] = starlet

            self.cluster_spectra_num[0] += len_spectra
            self.cluster_size_dist[0][len_spectra] = self.cluster_size_dist[0].get(len_spectra,0) + 1

        # build the noeds for cluster list1, and get the
        # distribution of cluster size.
        for cluster in self.cluster_lists[1]:
            len_spectra =  len(cluster.get_spectra())
            if len0 <= len1:
                starlet = ClusterNode(cluster.id, len_spectra) 
                self.starlets[cluster.id] = starlet
            else:
                star = ClusterNode(cluster.id, len_spectra) 
                self.stars[cluster.id] = star

            self.cluster_spectra_num[1] += len_spectra
            self.cluster_size_dist[1][len_spectra] = self.cluster_size_dist[1].get(len_spectra,0) + 1

        # do the comparing, and network building
        for i in range (0, len(self.cluster_lists[shorter_index])):
            cluster0 = self.cluster_lists[shorter_index][i] 
            for j in range (i, len(self.cluster_lists[longer_index])):
                cluster1 = self.cluster_lists[longer_index][j] 
                (shared_spec_num, similarity) = self.calculate_similarity(cluster0, cluster1)
                if similarity == 0:
                    continue
                self.similarity_dist[str(int(similarity*10))] = self.similarity_dist.get(str(int(similarity*10)),0) + 1
                self.shared_spec_num += shared_spec_num
                if len0 <= len1:
                    self.stars[cluster0.id].add_nb_node(cluster1.id, similarity, shared_spec_num)
                    self.starlets[cluster1.id].add_nb_node(cluster0.id, similarity, shared_spec_num)
                else:
                    self.stars[cluster1.id].add_nb_node(cluster0.id, similarity, shared_spec_num)
                    self.starlets[cluster0.id].add_nb_node(cluster1.id, similarity,shared_spec_num)

    
    def caculate_network_statistics(self):
        """
        Calculate the divide factors in the network.
        """
        self.standalone_star_num += 1 
        for key in self.stars.keys():
            star = self.stars[key]
            if star.nb_num == 0 :
                self.standalone_star_num += 1 
                continue
            divide_factor = star.nb_num + star.spec_num - star.shared_spec_num
            self.star_divide_factor_dist[divide_factor] = self.star_divide_factor_dist.get(divide_factor,0) + 1
            if star.spec_num < star.shared_spec_num:
                print("!!!!!!!!!!!!Becareful, total spectra No is less than Shared Spectra with starlets")
                print("with star " + star.id + "  " + str(star.spec_num) + "is less than" + str(star.shared_spec_num))
            if star.spec_num > star.shared_spec_num:
                star_lost_spec_numt += star.spec_num - star.shared_spec_num

        self.standalone_starlet_num += 1 
        for key in self.starlets.keys():
            starlet = self.starlets[key]
            if starlet.nb_num == 0 :
                self.standalone_starlet_num += 1 
                continue
            divide_factor = starlet.nb_num + starlet.spec_num - starlet.shared_spec_num
            self.starlet_divide_factor_dist[divide_factor] = self.starlet_divide_factor_dist.get(divide_factor,0) + 1
            if starlet.spec_num < starlet.shared_spec_num:
                print("!!!!!!!!!!!!Becareful, total spectra No is less than Shared Spectra with starlets")
                print("with star " + starlet.id + "  " + str(starlet.spec_num) + "is less than" + str(starlet.shared_spec_num))
            if starlet.spec_num > starlet.shared_spec_num:
                starlet_lost_spec_numt += starlet.spec_num - starlet.shared_spec_num


    def prepare_statistics(self):
        """
        Prepare the statistics output in row strings.
        """

        # statistics of clustering files
        percentage0 = "%.2f" % (100.0 * float(self.shared_spec_num)/float(self.cluster_spectra_num[0]))
        percentage1 = "%.2f" % (100.0 * float(self.shared_spec_num)/float(self.cluster_spectra_num[1]))
        head = "{0:<25}{1:<20}{2:<20}\n".format("name", "number", "description")
        rows = ""
        rows += "{0:<25}{1:<20}{2:<20}\n".format("cluster No.", len(self.cluster_lists[0]), "in file0")
        rows += "{0:<25}{1:<20}{2:<20}\n".format("cluster No.", len(self.cluster_lists[1]), "in file1")
        rows += "{0:<25}{1:<20}{2:<20}\n".format("identical cluster No.", self.similarity_dist['10'], "between them")
        rows += "{0:<25}{1:<20}{2:<20}\n".format("spectrum No", self.cluster_spectra_num[0], "in file0")
        rows += "{0:<25}{1:<20}{2:<20}\n".format("spectrum No", self.cluster_spectra_num[1], "in file1")
        rows += "{0:<25}{1:<20}{2:<20}\n".format("shared spectrum No", self.shared_spec_num, "between them")
        rows += "{0:<25}{1:<20}{2:<20}\n".format("shared spectrum percent", percentage0, "in file0")
        rows += "{0:<25}{1:<20}{2:<20}\n".format("shared spectrum percent", percentage1, "in file1")
        self.tables.append(('statistics of files', head, rows))

        # distribution of cluster size in file0
        head = '{0:<20}{1:<20}{2:<20}{3:<20}\n'.format("cluster size","No.", "percentage", "accumulate pecentage")
        rows = ""
        len0 = len(self.cluster_lists[0])
        accumulate_num = 0
        for key in self.cluster_size_dist[0].keys():
            value = self.cluster_size_dist[0][key]
            accumulate_num += value
            percent = "%.2f" % (100 * value/len0)
            accum_percent = "%.2f" % (100 * accumulate_num/len0)
            rows += '{0:<20}{1:<20}{2:<20}{3:<20}\n'.format(key, value, percent, accum_percent)
        self.tables.append(('distribution of cluster size in file0', head, rows))
        
        head = '{0:<20}{1:<20}{2:<20}{3:<20}\n'.format("cluster size","No.", "percentage", "accumulate pecentage")
        rows = ""
        len1 = len(self.cluster_lists[1])
        accumulate_num = 0
        for key in self.cluster_size_dist[1].keys():
            value = self.cluster_size_dist[1][key]
            accumulate_num += value
            percent = "%.2f" % (100 * value/len1)
            accum_percent = "%.2f" % (100 * accumulate_num/len1)
            rows += '{0:<20}{1:<20}{2:<20}{3:<20}\n'.format(key, value, percent, accum_percent)
        self.tables.append(('distribution of cluster size in file1', head, rows))

        # distribution of similarity
        head = "{0:<20}{1:<20}{2:<20}{3:<20}\n".format("similarity score", "pairs of clusters", "percentage(stars)", "percentage(starlets)")
        rows = ""
        for key in self.similarity_dist.keys():
            value = self.similarity_dist[key]
#            rows += "%d %s %d %s %.2f %s %.2f\n" % (int(key),  int(value),  100.0*value/self.stars_length,  100.0*value/self.starlets_length)
            percent_star = "%.2f" % (100.0*value/self.stars_length)
            percent_starlet  = "%.2f" % (100.0*value/self.starlets_length)
            rows += '{0:<20}{1:<20}{2:<20}{3:<20}\n'.format(key, value, percent_star, percent_starlet)
        self.tables.append(('distribution of similarity (identical = 10)', head, rows))

        # distribution of star divide factors
        head = '{0:<20}{1:<20}\n'.format("divide factor","No.")
        rows = ""
        for key in self.star_divide_factor_dist.keys():
            percent_star = "%.2f" % (100.0*value/self.stars_length)
            value = self.star_divide_factor_dist[key]
            rows += '{0:<20}{1:<20}{1:<20}\n'.format(key, value, percent_star)
        self.tables.append(('distribution of star divide factors', head, rows))

        # distribution of starlet divide factors
        head = '{0:<20}{1:<20}\n'.format("divide factor","No.")
        rows = ""
        for key in self.starlet_divide_factor_dist.keys():
            value = self.starlet_divide_factor_dist[key]
            percent_starlet  = "%.2f" % (100.0*value/self.starlets_length)
            rows += '{0:<20}{1:<20}{1:<20}\n'.format(key, value, percent_starlet)
        self.tables.append(('distribution of starlet divide factors', head, rows))


    
    def output_debug_info(self):
        """
        Only for debug, output the intersted intermediate data.
        """
        print(len(self.cluster_lists[0]))
        print(len(self.cluster_lists[1]))
#        print(self.tables['similarity_dist'])
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
