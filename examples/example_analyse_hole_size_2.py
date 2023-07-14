import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import networkx as nx

from capsidgraph.generator import (
    icosahedral_patterns,
    create_icosahedral_face_edges,
    create_icosahedral_capsid_graph,
)
from capsidgraph.analyser import (
    probability_fragment,
    energy_nodes_fragment,
    get_hole_size_distribution,
)

from capsidgraph.analyser import init_nodes_energy

iterations = 10000
frames = 30

h = 1
k = 1
[edges, Tx, Ty, Tscale] = icosahedral_patterns.PATTERN_666
face_edges, axis = create_icosahedral_face_edges(edges, Tx, Ty, h, k)
G = create_icosahedral_capsid_graph(face_edges, axis)
nx.set_edge_attributes(G, 1 / len(G.edges), "energy")
init_nodes_energy(G)
distributions = []
for p in np.linspace(0, 1, frames):
    print("Computing for p=", p)
    # distributions.append(get_hole_size_distribution(G, iterations, probability_fragment, {"fragmentation":p, "fragmentation_type":"nodes"}))
    distributions.append(
        get_hole_size_distribution(
            G, iterations, energy_nodes_fragment, {"fragmentation": p}
        )
    )

fig, ax = plt.subplots()

x = np.arange(0, len(G.nodes) + 1, 1)
bar = ax.bar(x, np.zeros(len(G.nodes) + 1))
ax.set_ylim(0, max([max(d) for d in distributions]))


def animate(i):
    ax.text(
        0.5,
        1.100,
        "Holes size distribution, p=%d" % (100 * i / frames),
        bbox={"facecolor": "white", "alpha": 1, "pad": 5},
        transform=ax.transAxes,
        ha="center",
    )
    for j in range(len(G.nodes) + 1):
        bar[j].set_height(distributions[i][j])
    return bar

plt.show()
