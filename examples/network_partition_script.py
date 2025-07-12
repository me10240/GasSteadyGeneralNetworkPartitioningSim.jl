import json
import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG) # decides which events should be propagated
import networkx as nx
# from networkx.algorithms.components import is_biconnected
from networkx.algorithms.community import kernighan_lin_bisection


class CustomLogFormatter(logging.Formatter):
    format = "%(asctime)s - %(name)s  - %(levelname)s - (%(filename)s:%(lineno)d  %(funcName)s) - %(message)s"
    FORMATS = {
        logging.DEBUG: format,
        logging.INFO: format,
        logging.WARNING: format,
        logging.ERROR: format,
        logging.CRITICAL: format
    }
    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt, "%Y-%m-%d %H:%M:%S")
        return formatter.format(record)


class CustomStreamFormatter(logging.Formatter):

    green = "\x1b[32;20m"
    grey = "\x1b[38;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format = "%(asctime)s - %(name)s - %(levelname)s - (%(filename)s:%(lineno)d  %(funcName)s) - %(message)s "

    FORMATS = {
        logging.DEBUG: green + format + reset,
        logging.INFO: grey + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt, "%Y-%m-%d %H:%M:%S")
        return formatter.format(record)


def load_data_and_create_graph(filename, plotting_flag=False):
    with open(filename, "r") as read_file:
        network_data = json.load(read_file)

    slack_nodes = []

    G = nx.Graph()
    G.clear()
    comp_list = ['pipes', 'short_pipes', 'compressors', 'valves', 'control_valves', 'resistors','loss_resistors']
    for comp in comp_list:
        if comp in network_data:
            for id in network_data[comp]:
                if "fr_node" in network_data[comp][id].keys():
                    G.add_edge(network_data[comp][id]['fr_node'], network_data[comp][id]['to_node'])
                else:
                    G.add_edge(network_data[comp][id]['from_node'], network_data[comp][id]['to_node'])

    for node_id in network_data['nodes']:
        G.nodes[int(node_id)]["pos"] = (network_data['nodes'][node_id]['x_coord'], network_data['nodes'][node_id]['y_coord'])
        if network_data['nodes'][node_id]['slack_bool'] == 1:
            slack_nodes.append(int(node_id))

    log.info("Network has {} nodes, {} edges, and the slack node(s) is/are {}".format(G.number_of_nodes(), G.number_of_edges(), slack_nodes))

    if plotting_flag:
        pos = nx.get_node_attributes(G, "pos")
        import matplotlib.pyplot as plt
        plt.figure(figsize=(4.5, 4.5), dpi=200)
        color_list = ["seagreen" if node_name in slack_nodes  else 'darkorange' for node_name in list(G.nodes)]
        nx.draw_networkx(G, pos=pos, with_labels=False, node_size= 5, font_size=5, font_weight='bold', node_color=color_list, edge_color="black")
        plt.show()
    return G, slack_nodes


def partition_graph_into_subgraphs_at_chosen_articulation_points(G, interface_nodes, allow_singletons_flag=False):
   
    if type(interface_nodes) == int:
        interface_node_list = [interface_nodes]
    else:
        interface_node_list = interface_nodes

    log.debug("interface_node_list:{}".format(interface_node_list))
    nbr_sets = []
    for i in interface_node_list:
        nbr_sets.append( list(G.neighbors(i)) )
    log.debug("nbr_sets:{}".format(nbr_sets))


    for i in range(len(interface_node_list)):
        nbrs = nbr_sets[i]
        node = interface_node_list[i]
        log.debug("node:{}, nbrs:{}".format(node, nbrs))
        for i in nbrs:
            if (node, i) in G.edges():
                G.remove_edge(node, i)

    len_subnetworks_array = [len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]

    
    for num in len_subnetworks_array:
        if num == 1:
            len_subnetworks_array.remove(num)

    if len(len_subnetworks_array) >= 2:
        log.info("Given node list {} successfully partitions network".format(interface_node_list))
    else:
        log.error("Given node list {} does not partition the network!".format(interface_node_list))
        exit()

    S = [G.subgraph(c).copy() for c in sorted(nx.connected_components(G), key=len, reverse=False)]

    removal_list = []
    for SG in S:
        for node in interface_node_list:
            if node in SG.nodes():
                removal_list.append(SG)
                break
    for SG in removal_list:
        S.remove(SG)

    if not allow_singletons_flag:
        # check if any singletons remain
        if S[0].number_of_nodes() == 1:
            log.error("Partitioning issue! Singletons found !")
            exit()

        
    for ni, SG in enumerate(S):
        for i, intf_node in enumerate(interface_node_list): 
            nbrs = nbr_sets[i]
            for node in nbrs:
                if node in SG.nodes():
                    SG.add_edge(intf_node, node) 
    return S

def partition_graph(G, interface_nodes, slack_nodes, allow_slack_node_partitioning=True):

    log.debug("Nodes:{}, edges:{}".format(G.number_of_nodes(), G.number_of_edges()))

    if not allow_slack_node_partitioning:
        interface_nodes = [node for node in interface_nodes if node not in slack_nodes] # to preserve sorted order
        if interface_nodes == []:
            log.warning("Given nodes are all slack nodes! Could not partition.")
            return [], []

    S  = partition_graph_into_subgraphs_at_chosen_articulation_points(G, interface_nodes, allow_singletons_flag=True)
    
    return S
    

def detect_slack_networks(S, slack_nodes):

    def slack_ordering(SG, slack_nodes):
        return len(set(SG.nodes()).intersection(set(slack_nodes)))
        
    S.sort(reverse=True, key=lambda SG: slack_ordering(SG, slack_nodes))

    slack_network_ids = []
    # identify slack networks
    for slack_id in slack_nodes:
        for ni, SG in enumerate(S):
            if slack_id in SG.nodes() and ni not in slack_network_ids:
                slack_network_ids.append(ni)
                log.debug("Found slack node {} in subnetwork {}+1 ".format(slack_id, ni))
    
    if len(slack_network_ids) == 1:
        log.info("All slack nodes within single subnetwork")
    else:
        log.warning("There are {} slack subnetworks".format(len(slack_network_ids)))

    return S, slack_network_ids

def write_partition_json_file_for_julia(filename, S, interface_node_list, slack_nodes):
    partition_dict = {}
    partition_dict["num_partitions"] = len(S)
    partition_dict["interface_nodes"] = interface_node_list
    partition_dict["slack_nodes"] = slack_nodes
    for ni, SG in enumerate(S):
        partition_dict[ni+1] = list(SG.nodes()) #starting from 1 instead of 0
    log.info("Writing to json...")
    with open(filename, "w") as outfile:
        json.dump(partition_dict, outfile)
    log.info("Completed writing to json...")
    return partition_dict

def construct_subnetwork_graph(p_dict, plotting_flag=False, full_network=None):

    G = nx.Graph()
    num_partitions = p_dict["num_partitions"]
    interface_node_list = p_dict["interface_nodes"]

    for i in range(1, num_partitions+1):
        for j in interface_node_list:
            if j in set(p_dict[i]):
                node = "N-{}".format(i)
                G.add_edge(node, j)
    tree_status = nx.is_tree(G)
    log.info("Block cut graph is a tree: {}".format(tree_status))

    if plotting_flag:
        if full_network:
            def compute_partition_centroid(G0, p_dict, i):
                centroid_x = 0.0
                centroid_y = 0.0
                num_nodes = len(p_dict[i])
                for node in p_dict[i]:
                    centroid_x += G0.nodes[node]["pos"][0]
                    centroid_y += G0.nodes[node]["pos"][1]
                return (centroid_x/num_nodes, centroid_y/num_nodes)
            
            for i in range(1, num_partitions+1):
                G.nodes[f"N-{i}"]["pos"] = compute_partition_centroid(full_network, p_dict, i)
            for i in interface_node_list:
                G.nodes[i]["pos"] = full_network.nodes[i]["pos"]
            pos = nx.get_node_attributes(G, "pos")
        else:
            pos = nx.spring_layout(G)

        # slack_network_names = [f"N-{i}" for i in slack_network_ids]
        import matplotlib.pyplot as plt
        plt.figure(figsize=(4.5, 4.5), dpi=200)
        color_list = ['tan' if node_name in interface_node_list   else 'darkorange' for node_name in list(G.nodes)]
        size_list = [50 if node_name in interface_node_list else 150 for node_name in list(G.nodes)]
        font_size_list = [2 if node_name in interface_node_list else 4 for node_name in list(G.nodes)]

        nx.draw_networkx(G, pos=pos, with_labels=True, node_size= size_list, font_size=5, font_weight='bold', node_color=color_list, edge_color="black", style="--")
        plt.show()

    return

def partition_given_network(dirname, node_list, allow_slack_node_partitioning=True, plotting_flag=False):
    filename = dirname + "network.json"
    G, slack_nodes = load_data_and_create_graph(filename, plotting_flag=plotting_flag)
    interface_node_list = node_list

    # check that there are no edges between nodes in interface_node_list
    def filter_node_fn(node):
        return node in interface_node_list
    
    Gsub = nx.subgraph_view(G, filter_node=filter_node_fn)
    if Gsub.number_of_edges() > 0:
        log.error("There are edges between nodes in interface_node_list")
        exit()
    else:
        log.info("No edges between nodes in {}".format(interface_node_list))


    S = partition_graph(G, interface_node_list, slack_nodes, allow_slack_node_partitioning=allow_slack_node_partitioning)

    S, slack_network_ids = detect_slack_networks(S, slack_nodes)

    partition_sizes = [SG.number_of_nodes()   for SG in S ]
    
    log.info("Size of largest partition is {}".format(max(partition_sizes)))

    if len(S) >= 2:
        partition_data_file = dirname + "partition-test-script-dummy.json"
        partition_dict = write_partition_json_file_for_julia(partition_data_file, S, interface_node_list, slack_nodes)
        construct_subnetwork_graph(partition_dict, plotting_flag=plotting_flag, full_network=None)
    else:
        log.info("Could not partition network")
    return

def main():
    
    import os
    print(os.getcwd())

    # dirname = "./data/Texas7k_Gas/"
    dirname = "./data/GasLib-40/"
    node_list = [17, 27, 32]

    run_script(dirname, node_list, loglevel="info", allow_slack_node_partitioning = False, plotting_flag=True)


def run_script(dirname, node_list, loglevel="info", allow_slack_node_partitioning = False, plotting_flag=True):
    
    
    level = getattr(logging, loglevel.upper())  #levels are 10, 20, 30, 40, 50

    # create console handler
    ch = logging.StreamHandler()
    ch.setLevel(level)
    ch.setFormatter(CustomStreamFormatter()) 
    log.addHandler(ch)
    
    # create file handler
    logfile = dirname + "partition_script.log"
    fh = logging.FileHandler(logfile, mode='w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(CustomLogFormatter())
    log.addHandler(fh) 

    log.info("Using NetworkX version {}".format(nx.__version__))
    # dirname = "8-node/"
    partition_given_network(dirname, node_list, allow_slack_node_partitioning = allow_slack_node_partitioning, plotting_flag=plotting_flag)

    
if __name__ == "__main__":
    main()