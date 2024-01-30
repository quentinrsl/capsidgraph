import networkx as nx
from capsidgraph.generator import icosahedral_patterns,create_icosahedral_face_edges,create_icosahedral_capsid_graph
from capsidgraph.analyser import *

# This is an example of how the bisection could be used to compute fragmentation threshold with a custom fragmentation function
# Here instead of considering a fragmented graph as a graph with two connected components, we consider a graph as fragmented if it has a hole of size at least half of the original graph
# The method to get the threshold is the same as in the other examples, we find the parameter that makes the probability of fragmentation equal to 1/2 with a bisection method

steps = 7
error_prob = 0.1
max_iterations = 1000000
process_number = 4

h = 2
k = 1
[edges, Tx, Ty, Tscale] = icosahedral_patterns.PATTERN_333333
face_edges, axis = create_icosahedral_face_edges(edges, Tx, Ty, h, k)
G_original = create_icosahedral_capsid_graph(face_edges, axis)
nx.set_edge_attributes(G_original, 1/len(G_original.edges), "strength")
init_nodes_strength(G_original)

def is_fragmented(G,G_original=G_original):
    s = get_hole_size(G,G_original)
    return s >= len(G_original.nodes)//2

hole_unweighted_threshold_nodes = bisection(G_original, steps, error_prob, probability_fragment,is_fragmented=is_fragmented, fragment_settings={"fragmentation_type" : "nodes"},debug=True, max_iterations=max_iterations, process_number=process_number)
hole_unweighted_threshold_edges = bisection(G_original, steps, error_prob, probability_fragment,is_fragmented=is_fragmented, fragment_settings={"fragmentation_type" : "edges"},debug=True, max_iterations=max_iterations, process_number=process_number)
hole_weighted_threshold_nodes = bisection(G_original, steps, error_prob, strength_nodes_fragment,is_fragmented=is_fragmented, fragment_settings={},debug=True, max_iterations=max_iterations, process_number=process_number)

print("Unweighted threshold nodes: ", hole_unweighted_threshold_nodes)
print("Unweighted threshold edges: ", hole_unweighted_threshold_edges)
print("Weighted threshold nodes: ", hole_weighted_threshold_nodes)