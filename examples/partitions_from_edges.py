import json
import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG) # decides which events should be propagated
import networkx as nx
from itertools import product

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

def load_data_and_create_graph(network_data):
    
    slack_nodes = []
    G = nx.Graph()
    G.clear()
    comp_list = ['pipes', 'short_pipes', 'compressors', 'valves', 'control_valves', 'resistors','loss_resistors']
    for comp in comp_list:
        if comp in network_data:
            for id in network_data[comp]:
                fr_str = "fr_node" if "fr_node" in network_data[comp][id].keys() else "from_node"
                to_str = "to_node"
                fr_node = network_data[comp][id][fr_str]
                to_node = network_data[comp][id][to_str]
                if G.has_edge(fr_node, to_node):
                    log.warning("MORE THAN ONE EDGE BETWEEN NODES {} AND {}".format(fr_node, to_node))
                    continue
                G.add_edge(fr_node, to_node)
                

    for node_id in network_data['nodes']:
        if network_data['nodes'][node_id]['slack_bool'] == 1:
            slack_nodes.append(int(node_id))

    log.info("Network has {} nodes, {} edges, and the slack node(s) is/are {}".format(G.number_of_nodes(), G.number_of_edges(), slack_nodes))

    return G, slack_nodes

def modified_network(G, edge_removal):
    remove_edge_list = []
    for _, edge in enumerate(edge_removal):
        to_node = edge["to_node"]
        fr_node = edge.get("from_node", edge.get("fr_node"))
        remove_edge_list.append((fr_node, to_node))
    edges = [(u,v) for (u,v) in list(G.edges) if (u,v) not in remove_edge_list and (v,u) not in remove_edge_list]
    G_truncated =  nx.edge_subgraph(G, edges)
    return G_truncated, remove_edge_list

def create_partition_data(G, G_truncated, remove_edge_list, slack_nodes):
    S = nx.connected_components(G_truncated)
    partitions = [list(c) for c in S]
    def slack_ordering(s, slack_nodes):
        return len(set(s).intersection(set(slack_nodes)))
        
    partitions.sort(reverse=True, key=lambda c: slack_ordering(list(c), slack_nodes))
    num_partitions = len(partitions)

    partition_dict = {}
    partition_dict["num_partitions"] = num_partitions
    partition_dict["slack_nodes"] = slack_nodes
    partition_dict["interface_nodes"] = []

    #construct vertex sequence
    chosen_seq = None
    for seq in product(*remove_edge_list): # * is the unpacking operator 
        num_edge = nx.number_of_edges(nx.induced_subgraph(G, set(seq)))
        if num_edge == 0:
            chosen_seq = set(seq)
            break
    if chosen_seq is None:
        log.error("Unable to find  a vertex separator corresponding to given edge separator")
        exit()


    for i in chosen_seq:
        partition_dict["interface_nodes"].append(i)
        for (u,v) in remove_edge_list:
            if i == u or i == v:
                vertex = v if i == u else u
                for k in range(len(partitions)):
                    if i in partitions[k] and v not in partitions[k]:
                        partitions[k].append(vertex)
    for i in range(len(partitions)):
        partition_dict[i+1] = partitions[i]

    slack_network_ids = []
    # identify slack networks
    for slack_id in slack_nodes:
        for ni, c in enumerate(partitions):
            if slack_id in list(c) and ni not in slack_network_ids:
                slack_network_ids.append(ni+1)
    partition_dict["slack_network_ids"] = slack_network_ids            
    return partition_dict
    
def write_partition_json_file(filename, partition_dict):
    
    log.info("Writing to json...")
    with open(filename, "w") as outfile:
        json.dump(partition_dict, outfile)
    log.info("Completed writing to json...")
    return

def run_script(dirname, partition_file, loglevel="info"):
    
    
    level = getattr(logging, loglevel.upper())  #levels are 10, 20, 30, 40, 50

    # create console handler
    ch = logging.StreamHandler()
    ch.setLevel(level)
    ch.setFormatter(CustomStreamFormatter()) 
    log.addHandler(ch)
    
    # create file handler
    logfile = dirname + "partitions_from_edges.log"
    fh = logging.FileHandler(logfile, mode='w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(CustomLogFormatter())
    log.addHandler(fh) 

    log.info("Using NetworkX version {}".format(nx.__version__))
    network_filename = dirname + "network.json"
    with open(network_filename, "r") as read_file:
        network_data = json.load(read_file)
    G, slack_nodes = load_data_and_create_graph(network_data)

    # piope[14] 17 --32
    # pipe[15] 36--17
    # pipe[39] 32--27
    # pipe[36] 26--27


    # pipe[17] 36--21
    # pipe[2] 32--21
    
    # edge_removal =[network_data['pipes']['3'], network_data['compressors']['3'],network_data['compressors']['1']] # for GasLib24
    edge_removal =[network_data['pipes']['14'], network_data['pipes']['15'], network_data['pipes']['39'], network_data['pipes']['36']] # for GasLib40

    G_truncated, remove_edge_list = modified_network(G, edge_removal)
    if not nx.is_connected(G_truncated):
        partition_dict = create_partition_data(G, G_truncated, remove_edge_list, slack_nodes)
        write_partition_json_file(dirname+partition_file, partition_dict)

def main():
    
    import os
    print(os.getcwd())

    dirname = "./data/GasLib-40/"
    partition_file = "partition_data_delete.json"

    run_script(dirname, partition_file, loglevel="info")
    

if __name__ == "__main__":
    main()