import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

from capsidgraph.generator import (
    icosahedral_patterns,
    create_icosahedral_face_edges,
    create_icosahedral_capsid_graph,
)

for P in [icosahedral_patterns.PATTERN_333333, icosahedral_patterns.PATTERN_666, icosahedral_patterns.PATTERN_6363, icosahedral_patterns.PATTERN_6434, icosahedral_patterns.PATTERN_3336]:
    h = 1
    k = 1
    [edges, Tx, Ty, Tscale] = P
    face_edges, axis = create_icosahedral_face_edges(edges, Tx, Ty, h, k)
    G = create_icosahedral_capsid_graph(face_edges, axis)

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
