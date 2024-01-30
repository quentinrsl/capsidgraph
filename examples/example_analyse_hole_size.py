import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

from capsidgraph.generator import (
    icosahedral_patterns,
    create_icosahedral_face_edges,
    create_icosahedral_capsid_graph,
)
from capsidgraph.analyser import (
    strength_nodes_fragment,
    init_nodes_strength,
    get_hole_size_distribution
)

"""
This example shows how to use the get_hole_size_distribution function 
to determine the probability distribution of hole sizes of a graph under fargmentation.
The fragmentation method used here is the strength node removal method with all edges having the same strength.
the result is diplayed in a matplotlib bar graph
"""

iterations = 10000
E = 1

h = 1
k = 1
[edges, Tx, Ty, Tscale] = icosahedral_patterns.PATTERN_666
face_edges, axis = create_icosahedral_face_edges(edges, Tx, Ty, h, k)
G = create_icosahedral_capsid_graph(face_edges, axis)
nx.set_edge_attributes(G, 1 / len(G.edges), "strength")
init_nodes_strength(G)
res = get_hole_size_distribution(G, iterations, strength_nodes_fragment, {"fragmentation":E, "fragmentation_type":"nodes"})

plt.bar(list(range(len(res))), res)
plt.show()