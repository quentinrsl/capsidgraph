import matplotlib.pyplot as plt

from capsidgraph.generator import (
    icosahedral_patterns,
    create_icosahedral_face_edges,
    create_icosahedral_capsid_graph,
)
from capsidgraph.analyser import probability_fragment, get_fragment_size_distribution

"""
This example shows how to use the get_fragment_size_distribution function 
to determine the average fragment size distribution of a graph under fargmentation.
This example uses random node removal and displays the result in a matplotlib bar graph.
"""

iterations = 10000
p = 0.2

h = 1
k = 1
[edges, Tx, Ty, Tscale] = icosahedral_patterns.PATTERN_666
face_edges, axis = create_icosahedral_face_edges(edges, Tx, Ty, h, k)
G = create_icosahedral_capsid_graph(face_edges, axis)
res = get_fragment_size_distribution(
    G,
    iterations,
    probability_fragment,
    {"fragmentation": p, "fragmentation_type": "nodes"},
)

plt.bar(list(range(len(res))), res)
plt.show()
