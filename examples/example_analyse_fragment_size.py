import matplotlib.pyplot as plt

from capsidgraph.generator import (
    icosahedral_patterns,
    create_icosahedral_face_edges,
    create_icosahedral_capsid_graph,
)
from capsidgraph.analyser import (
    probability_fragment,
    get_fragment_size_distribution
)

iterations = 10000

h = 1
k = 1
[edges, Tx, Ty, Tscale] = icosahedral_patterns.PATTERN_666
face_edges, axis = create_icosahedral_face_edges(edges, Tx, Ty, h, k)
G = create_icosahedral_capsid_graph(face_edges, axis)
res = get_fragment_size_distribution(G, iterations, probability_fragment, {"fragmentation":0.22, "fragmentation_type":"nodes"})

plt.bar(list(range(len(res))), res)
plt.show()