import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

from capsidgraph.generator import (
    cubic_patterns,
    create_cubic_face_edges,
    create_cubic_capsid_graph,
)

"""
This example shows how to use capsidgraph.generator to create a cubic capsid graph from a lattice parttern.
"""

for P in [cubic_patterns.AALS_24_PATTERN, cubic_patterns.AALS_48_PATTERN, cubic_patterns.AALS_60_PATTERN]:
    [edges, Tx, Ty, face_side_edge] = P
    face_edges, face_square_vertices = create_cubic_face_edges(edges, Tx, Ty, face_side_edge)
    G = create_cubic_capsid_graph(face_edges, face_square_vertices)

    pos = nx.kamada_kawai_layout(G, dim=3)

    node_xyz = np.array([pos[v] for v in sorted(G)])
    edge_xyz = np.array([(pos[u], pos[v]) for u, v in G.edges()])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(*node_xyz.T, s=100, ec="w")
    for vizedge in edge_xyz:
        ax.plot(*vizedge.T, color="tab:gray")

    ax.grid(False)
    for dim in (ax.xaxis, ax.yaxis, ax.zaxis):
        dim.set_ticks([])

    fig.tight_layout()
    plt.show()
