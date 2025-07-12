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


def load_data_and_create_graph(network_filename, partition_dict):
    with open(network_filename, "r") as read_file:
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

    
    pos = nx.get_node_attributes(G, "pos")
    import matplotlib.pyplot as plt
    plt.figure(figsize=(4.5, 4.5), dpi=200)
    color_list = ["seagreen" if node_name in slack_nodes  else 'darkorange' for node_name in list(G.nodes)]
    nx.draw_networkx(G, pos=pos, with_labels=False, node_size= 5, font_size=5, font_weight='bold', node_color=color_list, edge_color="black")
    

    num_partitions = partition_dict["num_partitions"]
    color_array = ["darkorange", "seagreen", "cornflowerblue", "cyan", "magenta", "peru", "red", "salmon", "darkviolet", "springgreen"]
    plt.figure(figsize=(4.5, 4.5), dpi=200)
    color_list2 = []
    for node_name in list(G.nodes):
        if node_name in partition_dict["interface_nodes"]:
            color_list2.append("black")
        else:
            for i in range(1, num_partitions+1):
                if node_name in partition_dict[str(i)]:
                    color_list2.append(color_array[i-1])
                    break

    node_size_list = [9 if node_name in partition_dict["interface_nodes"]  else 6 for node_name in list(G.nodes)]               
    nx.draw_networkx_nodes(G, pos=pos, node_size = node_size_list, node_color=color_list2, alpha=0.5)
    nx.draw_networkx_edges(G, pos=pos, node_size = node_size_list, edge_color="gray", width=0.5, alpha=0.5)
    plt.show()
    return G, slack_nodes



def construct_subnetwork_graph(p_dict, full_network=None):

    G = nx.Graph()
    num_partitions = p_dict["num_partitions"]
    interface_node_list = p_dict["interface_nodes"]

    for i in range(1, num_partitions+1):
        for j in interface_node_list:
            if j in set(p_dict[str(i)]):
                node = "N-{}".format(i)
                G.add_edge(node, j)
    tree_status = nx.is_tree(G)
    log.info("Subnetwork  hypergraph is a tree: {}".format(tree_status))

    
    if full_network:
        def compute_partition_centroid(G0, p_dict, i):
            centroid_x = 0.0
            centroid_y = 0.0
            num_nodes = len(p_dict[str(i)])
            for node in p_dict[str(i)]:
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

    import matplotlib.pyplot as plt
    plt.figure(figsize=(4.5, 4.5), dpi=200)
    color_list = ['tan' if node_name in interface_node_list else 'darkorange' for node_name in list(G.nodes)]
    size_list = [50 if node_name in interface_node_list else 150 for node_name in list(G.nodes)]
    font_size_list = [2 if node_name in interface_node_list else 4 for node_name in list(G.nodes)]

    nx.draw_networkx(G, pos=pos, with_labels=True, node_size= size_list, font_size=5, font_weight='bold', node_color=color_list, edge_color="black", style= "--")
    plt.show()

    return

def partition_given_network(dirname, partition_file):
    
    network_filename = dirname + "network.json"
    partition_filename = dirname + partition_file

    with open(partition_filename, "r") as read_file:
        partition_dict = json.load(read_file)
        
    G, slack_nodes = load_data_and_create_graph(network_filename, partition_dict)

    interface_node_list = partition_dict["interface_nodes"]

    # check that there are no edges between nodes in interface_node_list
    def filter_node_fn(node):
        return node in interface_node_list
    
    Gsub = nx.subgraph_view(G, filter_node=filter_node_fn)
    if Gsub.number_of_edges() > 0:
        log.error("There are edges between nodes in interface_node_list")
        exit()
    else:
        log.info("No edges between nodes in {}".format(interface_node_list))


    num_partitions = partition_dict["num_partitions"]
    partition_sizes = [len(partition_dict[str(i)])   for i in range(1, num_partitions+1)]
    
    log.info("Size of largest partition is {}".format(max(partition_sizes)))

    
    construct_subnetwork_graph(partition_dict, full_network=G)
    
    return

def main():
    
    import os
    print(os.getcwd())

    dirname = "./data/Texas7k_Gas/"
    # dirname = "./data/GasLib-40/"
    partition_file = "partition_data.json"

    run_script(dirname, partition_file, loglevel="info")


def run_script(dirname, partition_file, loglevel="info"):
    
    
    level = getattr(logging, loglevel.upper())  #levels are 10, 20, 30, 40, 50

    # create console handler
    ch = logging.StreamHandler()
    ch.setLevel(level)
    ch.setFormatter(CustomStreamFormatter()) 
    log.addHandler(ch)
    
    # create file handler
    logfile = dirname + "visualize_partition.log"
    fh = logging.FileHandler(logfile, mode='w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(CustomLogFormatter())
    log.addHandler(fh) 

    log.info("Using NetworkX version {}".format(nx.__version__))
    # dirname = "8-node/"
    partition_given_network(dirname, partition_file)

    
if __name__ == "__main__":
    main()