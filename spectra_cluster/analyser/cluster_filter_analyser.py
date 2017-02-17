"""
"""
import re,sys
from . import common


class ClusterFilter(common.AbstractAnalyser):
    """
    """
    def __init__(self, result_file):
        """
        """
        super().__init__()

        self.result_file = result_file

    def process_cluster(self, cluster):
        """
        """
        if self._ignore_cluster(cluster):
            return

        lines = list()    
        lines.append("\n=Cluster=")
        lines.append("id=" + cluster.id)
        lines.append("av_precursor_mz=" + str(cluster.precursor_mz))
        lines.append("av_precursor_intens=1.0")
        line = ""
        for mz in cluster.consensus_mz:
            line += "," + str(mz)
        lines.append("consensus_mz=" + line[1:])
        
        line = ""
        for intens in cluster.consensus_intens:
            line += "," + str(intens) 
        lines.append("consensus_intens=" + line[1:])
        
        for spectrum in cluster.get_spectra():
            p = re.compile("^.*?(P\wD\d{6}\;.*\;.*(spectrum|index|scan)=\d+)$")      
            m = p.match(spectrum.title)
            if m:
                line = "SPEC" + "\t"  +  m.group(1) + "\t" + str(spectrum.precursor_mz) \
                    + "\t" + str(spectrum.charge)  + "\t" + ",".join(spectrum.taxids) 
            else:
                print("ERROR, can not retreive spectrum id from title:" + spectrum.title)
                sys.exit(1)

            lines.append(line)
            
        # write the cluster 
        self.result_file.write("\n".join(lines))

    def remove_spectra_by_project(self, cluster, project_ids):
        temp_spectra = list()
        for spectrum in cluster.get_spectra():
            spectrum_title = spectrum.get_title()           
            p = re.compile("^.*?(P\wD\d{6})\;.*\;.*(spectrum|index|scan)=\d+$")      
            m = p.match(spectrum.title)
            if m:
                project_acc = m.group(1) 
            else:
                print("ERROR, can not retreive spectrum id from title:" + spectrum.title)
                sys.exit(1)

            if project_acc in project_ids:
                temp_spectra.append(spectrum)
        cluster.set_spectra(temp_spectra)


