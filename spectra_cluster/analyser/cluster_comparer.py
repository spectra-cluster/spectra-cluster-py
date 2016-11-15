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
        self.star_lost_spec_num = 0
        self.starlet_lost_spec_num = 0       
        self.ave_divide_factor_star = 0
        self.ave_divide_factor_starlet = 0
        self.ave_star_size = 0
        self.ave_starlet_size = 0
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
            return (spec.get_title())

        spectra.sort(key = mixed_order) 
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
            title0 = spectrum0.get_title() 
            title1 = spectrum1.get_title() 
            if(title0 < title1): 
                return -1
            elif(title0 > title1):   
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
        return (n,similarity_score)


    def compare(self):
        """
        Compare two cluster lists, get the distribution of the similarities.
        And build the connection network between clusters.
        """
        len0 = len(self.cluster_lists[0])
        len1 = len(self.cluster_lists[1])
        longer_index = 0 if len0 >= len1 else 1
        shorter_index = 1 if len1 <= len0 else 0
        self.stars_length = len(self.cluster_lists[shorter_index]) 
        self.starlets_length = len(self.cluster_lists[longer_index]) 
        # build the noeds for shorter cluster list, and get the
        # distribution of cluster size.
        for cluster in self.cluster_lists[shorter_index]:
            len_spectra =  len(cluster.get_spectra())
            star = ClusterNode(cluster.id, len_spectra) 
            self.stars[cluster.id] = star

            self.cluster_spectra_num[shorter_index] += len_spectra
            self.cluster_size_dist[shorter_index][len_spectra] = self.cluster_size_dist[shorter_index].get(len_spectra,0) + 1
        # build the noeds for longer cluster list, and get the
        # distribution of cluster size.
        for cluster in self.cluster_lists[longer_index]:
            len_spectra =  len(cluster.get_spectra())
            starlet = ClusterNode(cluster.id, len_spectra) 
            self.starlets[cluster.id] = starlet

            self.cluster_spectra_num[longer_index] += len_spectra
            self.cluster_size_dist[longer_index][len_spectra] = self.cluster_size_dist[longer_index].get(len_spectra,0) + 1
        # do the comparing, and network building
        for i in range (0, len(self.cluster_lists[shorter_index])):
            cluster0 = self.cluster_lists[shorter_index][i] 
            for j in range (i, len(self.cluster_lists[longer_index])):
                cluster1 = self.cluster_lists[longer_index][j] 
                (shared_spec_num, similarity) = self.calculate_similarity(cluster0, cluster1)
                if similarity == 0:
                    continue
                self.similarity_dist[int(similarity*10)] = self.similarity_dist.get(int(similarity*10),0) + 1
                self.shared_spec_num += shared_spec_num

                self.stars[cluster0.id].add_nb_node(cluster1.id, similarity, shared_spec_num)
                self.starlets[cluster1.id].add_nb_node(cluster0.id, similarity, shared_spec_num)

        self.ave_star_size = self.cluster_spectra_num[shorter_index]/self.stars_length
        self.ave_starlet_size = self.cluster_spectra_num[longer_index]/self.starlets_length
    
    def caculate_network_statistics(self):
        """
        Calculate the divide factors in the network.
        """
        divide_factor_sum = 0 
        for key in self.stars.keys():
            star = self.stars[key]
            if star.nb_num == 0 :
                self.standalone_star_num += 1 

            divide_factor = star.nb_num + 2 * (star.spec_num - star.shared_spec_num )/self.ave_starlet_size
            divide_factor_sum += divide_factor
            divide_factor_int = round(divide_factor)
            self.star_divide_factor_dist[divide_factor_int] = self.star_divide_factor_dist.get(divide_factor_int,0) + 1
            if star.spec_num < star.shared_spec_num:
                print("!!!!!!!!!!!!Becareful, total spectra No is less than Shared Spectra with starlets")
                print("with star " + star.id + "  " + str(star.spec_num) + "is less than" + str(star.shared_spec_num))
            if star.spec_num > star.shared_spec_num:
                self.star_lost_spec_num += star.spec_num - star.shared_spec_num
        self.ave_divide_factor_star = divide_factor_sum/self.stars_length

        divide_factor_sum = 0 
        for key in self.starlets.keys():
            starlet = self.starlets[key]
            if starlet.nb_num == 0 :
                self.standalone_starlet_num += 1 

            divide_factor = starlet.nb_num + 2 * (starlet.spec_num - starlet.shared_spec_num )/self.ave_star_size
            divide_factor_sum += divide_factor
            divide_factor_int = round(divide_factor)
            self.starlet_divide_factor_dist[divide_factor_int] = self.starlet_divide_factor_dist.get(divide_factor_int,0) + 1
            if starlet.spec_num < starlet.shared_spec_num:
                print("!!!!!!!!!!!!Becareful, total spectra No is less than Shared Spectra with starlets")
                print("with star " + starlet.id + "  " + str(starlet.spec_num) + "is less than" + str(starlet.shared_spec_num))
            if starlet.spec_num > starlet.shared_spec_num:
                self.starlet_lost_spec_num += starlet.spec_num - starlet.shared_spec_num
        self.ave_divide_factor_starlet = divide_factor_sum/self.starlets_length

    def prepare_statistics(self):
        """
        Prepare the statistics output in row strings.
        """

        # statistics of clustering files
        len0 = len(self.cluster_lists[0])
        len1 = len(self.cluster_lists[1])
        longer_index = 0 if len0 >= len1 else 1
        shorter_index = 1 if len1 <= len0 else 0

        percentage_stars = "%.2f" % (100.0 * float(self.shared_spec_num)/float(self.cluster_spectra_num[shorter_index]))
        percentage_starlets = "%.2f" % (100.0 * float(self.shared_spec_num)/float(self.cluster_spectra_num[longer_index]))

        head = "{0:<25}{1:<20}{2:<20}\n".format("name", "number", "description")
        rows = ""
        rows += "{0:<25}{1:<20}{2:<20}\n".format("stars No.", self.stars_length, "in file with less(or equal) clusters: file" + str(shorter_index))
        rows += "{0:<25}{1:<20}{2:<20}\n".format("starlets No.", self.starlets_length, "in file with more(or equal) clusters: file" + str(longer_index))
        rows += "{0:<25}{1:<20}{2:<20}\n".format("identical cluster No.", self.similarity_dist[10], "between them")
        rows += "{0:<25}{1:<20}{2:<20}\n".format("spectrum No", self.cluster_spectra_num[shorter_index], "in stars")
        rows += "{0:<25}{1:<20}{2:<20}\n".format("spectrum No", self.cluster_spectra_num[longer_index], "in starlets ")
        rows += "{0:<25}{1:<20}{2:<20}\n".format("shared spectrum No", self.shared_spec_num, "between them")
        rows += "{0:<25}{1:<20}{2:<20}\n".format("shared spectrum percent", percentage_stars, "in stars")
        rows += "{0:<25}{1:<20}{2:<20}\n".format("shared spectrum percent", percentage_starlets, "in starlets")
        self.tables.append(('statistics of files', head, rows))

        # distribution of cluster size in stars
        head = '{0:<20}{1:<20}{2:<20}{3:<20}\n'.format("cluster size","No.", "percentage", "accumulate pecentage")
        rows = ""
        rows += "{0:<20}{1:<20}\n".format("%.2f" % (self.ave_star_size), "average")
        accumulate_num = 0
        for key in sorted(self.cluster_size_dist[shorter_index].keys()):
            value = self.cluster_size_dist[shorter_index][key]
            accumulate_num += value
            percent = "%.2f" % (100 * value/self.stars_length)
            accum_percent = "%.2f" % (100 * accumulate_num/self.stars_length)
            rows += '{0:<20}{1:<20}{2:<20}{3:<20}\n'.format(key, value, percent, accum_percent)
        self.tables.append(('distribution of cluster size in stars', head, rows))
        
        head = '{0:<20}{1:<20}{2:<20}{3:<20}\n'.format("cluster size","No.", "percentage", "accumulate pecentage")
        rows = ""
        rows += "{0:<20}{1:<20}\n".format("%.2f" % (self.ave_starlet_size), "average")
        accumulate_num = 0
        for key in sorted(self.cluster_size_dist[longer_index].keys()):
            value = self.cluster_size_dist[longer_index][key]
            accumulate_num += value
            percent = "%.2f" % (100 * value/self.starlets_length)
            accum_percent = "%.2f" % (100 * accumulate_num/self.starlets_length)
            rows += '{0:<20}{1:<20}{2:<20}{3:<20}\n'.format(key, value, percent, accum_percent)
        self.tables.append(('distribution of cluster size in starlets', head, rows))

        # distribution of similarity
        head = "{0:<20}{1:<20}{2:<20}{3:<20}\n".format("similarity score", "pairs of clusters", "percentage(stars)", "percentage(starlets)")
        rows = ""
        for key in reversed(sorted(self.similarity_dist.keys())):
            value = self.similarity_dist[key]
            percent_star = "%.2f" % (100.0*value/self.stars_length)
            percent_starlet  = "%.2f" % (100.0*value/self.starlets_length)
            rows += '{0:<20}{1:<20}{2:<20}{3:<20}\n'.format(key, value, percent_star, percent_starlet)
        self.tables.append(('distribution of similarity (identical = 10)', head, rows))

        # distribution of star divide factors
        head = '{0:<20}{1:<20}{2:<20}\n'.format("divide factor","No.","percentage")
        rows = ""
        rows += "{0:<20}{1:<20}\n".format("%.2f" % (self.ave_divide_factor_star), "average")
        for key in sorted(self.star_divide_factor_dist.keys()):
            value = self.star_divide_factor_dist[key]
            percent_star = "%.2f" % (100.0*value/self.stars_length)
            rows += '{0:<20}{1:<20}{2:<20}\n'.format(key, value, percent_star)
        self.tables.append(('distribution of star divide factors', head, rows))

        # distribution of starlet divide factors
        head = '{0:<20}{1:<20}{2:<20}\n'.format("divide factor","No.","percentage")
        rows = ""
        rows += "{0:<20}{1:<20}\n".format("%.2f" % (self.ave_divide_factor_starlet), "average")
        for key in sorted(self.starlet_divide_factor_dist.keys()):
            value = self.starlet_divide_factor_dist[key]
            percent_starlet  = "%.2f" % (100.0*value/self.starlets_length)
            rows += '{0:<20}{1:<20}{2:<20}\n'.format(key, value, percent_starlet)
        self.tables.append(('distribution of starlet divide factors', head, rows))


    
    def output_debug_info(self):
        """
        Only for debug, output the intersted intermediate data.
        """
#        print self.sorted_spectra_dict[]

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
