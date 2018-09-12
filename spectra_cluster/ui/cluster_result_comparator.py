"""
This tools compares two clustering results with each other.

The core comparison is done by creating a network representation of the clustering results. Every cluster in
each of the files is represented as a node. Spectra present in both clusters create the edges.

The network is rendered in Cytoscape. For this to work Cytoscape must be opened before launching the script.
Additionally, the resulting network is stored as a GraphML file.

Usage:
  cluster_result_comparator.py --output <network.xml> <result1.clustering> <result2.clustering>

Options:
  -o, --output=<network.xml>           Path to the network file (GraphML format) to be created.
  -h, --help                           Print this help message.
  -v, --version                        Print the current version.
"""

import os
import sys

import networkx as nx
from docopt import docopt
from py2cytoscape.data.cyrest_client import CyRestClient
from py2cytoscape.data.style import StyleUtil

# make the spectra_cluster packages available
sys.path.insert(0, os.path.abspath('..') + os.path.sep + "..")

import spectra_cluster.clustering_parser as clustering_parser
from spectra_cluster.objects import Cluster


def load_spectra_to_cluster(result_file: str, before_cluster_id : str = None):
    """
    Creates a dict holding the spectra ids as keys and the cluster ids the spectrum belongs to as
    value.
    :param result_file: Path to the .clustering file to process
    :param before_cluster_id: If set this string is prepended to every cluster id.
    :return: A dict with the spectrum id as key and the cluster id as value
    """
    parser = clustering_parser.ClusteringParser(result_file)

    spec_to_cluster = dict()

    for cluster in parser:
        for spectrum in cluster.get_spectra():
            if before_cluster_id is not None:
                spec_to_cluster[spectrum.get_id()] = before_cluster_id + cluster.id
            else:
                spec_to_cluster[spectrum.get_id()] = cluster.id

    return spec_to_cluster


def build_network(spectra_to_cluster_1: dict, spectra_to_cluster_2: dict, source1: str, source2: str):
    """
    Builds a network with cluster's as nodes joined together by spectra
    shared by both clusters
    :param spectra_to_cluster_1: Spectrum id -> cluster id dict
    :param spectra_to_cluster_2: Spectrum id -> cluster id dict
    :param source1: Name of source 1
    :param source2: Name of source 2
    :return: The networkx object representing this network
    """
    graph = nx.Graph()

    # add all clusters as nodes
    graph.add_nodes_from(set(spectra_to_cluster_1.values()), source="source1", source_label=source1)
    graph.add_nodes_from(set(spectra_to_cluster_2.values()), source="source2", source_label=source2)

    # add the edges
    for spec_id_1 in spectra_to_cluster_1.keys():
        if spec_id_1 in spectra_to_cluster_2:
            cluster_1 = spectra_to_cluster_1[spec_id_1]
            cluster_2 = spectra_to_cluster_2[spec_id_1]

            try:
                graph.edges[cluster_1, cluster_2]["weight"] += 1
            except KeyError:
                graph.add_edge(cluster_1, cluster_2, weight=1)

    return graph


def remove_identical_clusters(graph: nx.Graph):
    """
    Removes all nodes from the graph that are only connected to a single
    other node (ie. identical clusters)
    :param graph: The graph to manipulate
    :return: The number of removed nodes
    """
    identical_clusters = list()
    for node in graph.nodes:
        neighbors = list(graph.neighbors(node))

        if len(neighbors) != 1:
            continue

        # make sure the neighbor also only has 1 neighbor
        neighbors_neighbors = list(graph.neighbors(neighbors[0]))

        if len(neighbors_neighbors) == 1:
            identical_clusters.append(node)

    graph.remove_nodes_from(identical_clusters)

    return len(identical_clusters)


def add_node_properties(graph: nx.Graph, clustering_file: str, label: str):
    """
    Add cluster properties to the network.
    :param graph: The graph to adapt
    :param clustering_file: The clustering file to read the clusters from
    :param label: The label used to generate the node ids in the graph
    """
    reader = clustering_parser.ClusteringParser(clustering_file)

    for cluster in reader:
        try:
            cluster_id = label + cluster.id
            graph.nodes[cluster_id]["size"] = cluster.n_spectra
            graph.nodes[cluster_id]["id"] = cluster.id
            graph.nodes[cluster_id]["precursor"] = cluster.precursor_mz
            graph.nodes[cluster_id]["min_mz"] = min([s.precursor_mz for s in cluster.get_spectra()])
            graph.nodes[cluster_id]["max_mz"] = max([s.precursor_mz for s in cluster.get_spectra()])

            # sequence ratios as string
            sequence_count_string = list()
            sequence_counts = Cluster.calculate_sequence_counts(cluster.get_spectra(), True)
            for sequence in sequence_counts:
                sequence_count_string.append(sequence + ": " + str(sequence_counts[sequence]))

            graph.nodes[cluster_id]["sequences"] = ", ".join(sequence_count_string)

        except KeyError:
            pass


def display_in_cytoscape(graph, network_name: str):
    """
    Displays the graph in Cytoscape - if it is running
    :param graph: The graph to display
    :param network_name: Name of the network in Cytoscape
    """
    try:
        cy = CyRestClient()
        cytoscape_network = cy.network.create_from_networkx(network=graph, name=network_name, collection="Cluster comparison")

        # add the layout
        cy.layout.apply(name='force-directed', network=cytoscape_network)

        # add the style
        # Step 6: Create Visual Style as code (or by hand if you prefer)
        comparison_style = cy.style.create('Cluster Comparison')

        basic_settings = {
            # You can set default values as key-value pairs.
            'NODE_FILL_COLOR': '#6AACB8',
            'NODE_SIZE': 55,
            'NODE_BORDER_WIDTH': 0,
            'NODE_LABEL_COLOR': '#555555',

            'EDGE_WIDTH': 2,
            'EDGE_TRANSPARENCY': 100,
            'EDGE_STROKE_UNSELECTED_PAINT': '#333333',

            'NETWORK_BACKGROUND_PAINT': '#FFFFFF'
        }

        comparison_style.update_defaults(basic_settings)

        # Create some mappings
        min_cluster_size = min([graph.node[c]["size"] for c in graph.nodes])
        max_cluster_size = max([graph.node[c]["size"] for c in graph.nodes])
        cluster_size_to_node_size = StyleUtil.create_slope(min=min_cluster_size, max=max_cluster_size, values=(30, 100))
        comparison_style.create_continuous_mapping(column='size', vp='NODE_SIZE', col_type='Double',
                                                   points=cluster_size_to_node_size)

        comparison_style.create_discrete_mapping(column='source', vp="NODE_FILL_COLOR", col_type="String",
                                                 mappings={"source1": "#FC8D62", "source2": "#A6D854"})
        comparison_style.create_passthrough_mapping(column='weight', vp="EDGE_WIDTH", col_type="Integer")

        # apply the style
        cy.style.apply(style=comparison_style, network=cytoscape_network)

    except Exception as e:
        print("Please open Cytoscape to display the graph: " + str(e))


def get_number_of_spectra(graph):
    """
    Get the total number of spectra (property size) of the graph
    :param graph: The graph to process
    :return: The total number of spectra
    """
    return sum([graph.node[c]["size"] for c in graph.nodes])


def main():
    arguments = docopt(__doc__, version='cluster_result_comparator.py 0.1')

    output_name = arguments["--output"]
    result_file_1 = arguments["<result1.clustering>"]
    result_file_2 = arguments["<result2.clustering>"]

    # sanity check
    if os.path.isfile(output_name):
        print("Error: Result file exists")
        sys.exit(1)

    if not os.path.isfile(result_file_1):
        print("Error: Failed to find " + result_file_1)
        sys.exit(1)
    if not os.path.isfile(result_file_2):
        print("Error: Failed to find " + result_file_2)
        sys.exit(1)

    label1 = os.path.basename(result_file_1).replace(".clustering", "_")
    label2 = os.path.basename(result_file_2).replace(".clustering", "_")

    # load the spectra to cluster information
    print("Loading data from " + result_file_1 + "...")
    spectra_to_cluster_1 = load_spectra_to_cluster(result_file_1, before_cluster_id=label1)
    print("  " + str(len(spectra_to_cluster_1)) + " spectra loaded")

    print("Loading data from " + result_file_2 + "...")
    spectra_to_cluster_2 = load_spectra_to_cluster(result_file_2, before_cluster_id=label2)
    print("  " + str(len(spectra_to_cluster_2)) + " spectra loaded")

    # build the network
    network = build_network(spectra_to_cluster_1, spectra_to_cluster_2, source1=os.path.basename(result_file_1),
                            source2=os.path.basename(result_file_2))

    print("Created network with " + str(len(network.nodes)) + " nodes and " + str(len(network.edges)) + " edges")

    # add node properties
    print("Adding node properties...")
    add_node_properties(network, result_file_1, label1)
    add_node_properties(network, result_file_2, label2)

    # count the total number of spectra
    total_spectra = get_number_of_spectra(network)

    # remove all nodes that are only connected to a single other node (ie. identical clusters)
    removed_nodes = remove_identical_clusters(network)

    print("Removed " + str(removed_nodes) + " identical nodes from the graph. " + str(len(network.nodes)) +
          " nodes remaining.")

    remaining_spectra = get_number_of_spectra(network)
    print("  Affecting " + str(remaining_spectra) + "/" + str(total_spectra) + " (" +
          str(round(remaining_spectra / total_spectra * 100, ndigits=2)) + "%) spectra")

    # save the GraphML file
    nx.write_graphml_xml(network, output_name)
    print("Network (GraphML format) written to " + output_name)

    print("Displaying comparison network in Cytoscape (must be open)...")
    display_in_cytoscape(network, network_name=label1 + " vs. " + label2)


if __name__ == "__main__":
    main()
