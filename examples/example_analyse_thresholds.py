import networkx as nx
from capsidgraph.generator import icosahedral_patterns,create_icosahedral_face_edges,create_icosahedral_capsid_graph
from capsidgraph.analyser import *

# This example shows how to use the library to compute fragmentation thresholds for different disassembly models, here wieghted and unweighted node removal as well as unweighted edge removal. 

steps = 7
error_prob = 0.1
max_iterations = 1000000
process_number = 4

h = 2
k = 1
[edges, Tx, Ty, Tscale] = icosahedral_patterns.PATTERN_333333
face_edges, axis = create_icosahedral_face_edges(edges, Tx, Ty, h, k)
G = create_icosahedral_capsid_graph(face_edges, axis)
nx.set_edge_attributes(G, 1/len(G.edges), "energy")
init_nodes_energy(G)

fragmentation_unweighted_threshold_edges = get_fragmentation_probability_threshold_edge(G, error_prob, steps, max_iterations=max_iterations, process_number=process_number, debug=True)
fragmentation_unweighted_threshold_nodes = get_fragmentation_probability_threshold_node(G, error_prob, steps, max_iterations=max_iterations, process_number=process_number, debug=True)
fragmentation_weighted_threshold_nodes = get_fragmentation_energy_threshold_node(G, error_prob, steps, max_iterations=max_iterations, process_number=process_number, debug=True)

print("Fragmentation probability threshold (unweighted) for edges: ", fragmentation_unweighted_threshold_edges)
print("Fragmentation probability threshold (unweighted) for nodes: ", fragmentation_unweighted_threshold_nodes)
print("Fragmentation energy threshold (weighted) for nodes: ", fragmentation_weighted_threshold_nodes)
