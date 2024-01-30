import unittest
import networkx as nx
import numpy as np

from capsidgraph.generator import create_icosahedral_face_edges
from capsidgraph.generator import create_icosahedral_capsid_graph
from capsidgraph.generator import create_cubic_capsid_graph
from capsidgraph.generator import create_icosahedral_texture
from capsidgraph.generator.face.patterns import icosahedral_patterns

from capsidgraph.generator import (
    create_cubic_face_edges,
    cubic_patterns,
    create_cubic_capsid_graph,
)

from math import sqrt


def strength_edge_match(e1, e2):
    return e1["strength"] == e2["strength"]


class TestGeneration(unittest.TestCase):
    def test_create_icosahedral_face_edges(self):
        pattern = icosahedral_patterns.PATTERN_333333
        h = 1
        k = 2
        [edges, Tx, Ty, Tscale] = pattern
        face_edges, axis = create_icosahedral_face_edges(edges, Tx, Ty, h, k)
        self.assertEqual(
            face_edges,
            [
                ((1, 0), (1, 1)),
                ((1, 1), (1, 2)),
                ((1, -1), (1, 0)),
                ((2, 0), (2, 1)),
                ((2, -1), (2, 0)),
                ((1, 0), (0, 1)),
                ((1, 1), (0, 2)),
                ((2, 0), (1, 1)),
                ((2, -1), (1, 0)),
                ((3, -1), (2, 0)),
                ((1, 0), (0, 0)),
                ((1, 1), (0, 1)),
                ((2, 0), (1, 0)),
                ((2, 1), (1, 1)),
                ((3, 0), (2, 0)),
                ((1, 0), (1, -1)),
                ((1, 1), (1, 0)),
                ((1, 2), (1, 1)),
                ((2, 0), (2, -1)),
                ((2, 1), (2, 0)),
                ((0, 1), (1, 0)),
                ((0, 2), (1, 1)),
                ((1, 0), (2, -1)),
                ((1, 1), (2, 0)),
                ((2, 0), (3, -1)),
                ((0, 0), (1, 0)),
                ((0, 1), (1, 1)),
                ((1, 0), (2, 0)),
                ((1, 1), (2, 1)),
                ((2, 0), (3, 0)),
            ],
        )
        self.assertEqual(axis, ((0, 0), (1, 2), (3, -1)))

    def test_icosahedral_graph_generation(self):
        face_edges = [
            ((1, 0), (1, 1)),
            ((1, 1), (1, 2)),
            ((1, -1), (1, 0)),
            ((2, 0), (2, 1)),
            ((2, -1), (2, 0)),
            ((1, 0), (0, 1)),
            ((1, 1), (0, 2)),
            ((2, 0), (1, 1)),
            ((2, -1), (1, 0)),
            ((3, -1), (2, 0)),
            ((1, 0), (0, 0)),
            ((1, 1), (0, 1)),
            ((2, 0), (1, 0)),
            ((2, 1), (1, 1)),
            ((3, 0), (2, 0)),
            ((1, 0), (1, -1)),
            ((1, 1), (1, 0)),
            ((1, 2), (1, 1)),
            ((2, 0), (2, -1)),
            ((2, 1), (2, 0)),
            ((0, 1), (1, 0)),
            ((0, 2), (1, 1)),
            ((1, 0), (2, -1)),
            ((1, 1), (2, 0)),
            ((2, 0), (3, -1)),
            ((0, 0), (1, 0)),
            ((0, 1), (1, 1)),
            ((1, 0), (2, 0)),
            ((1, 1), (2, 1)),
            ((2, 0), (3, 0)),
        ]
        axis = ((0, 0), (1, 2), (3, -1))
        G1 = create_icosahedral_capsid_graph(face_edges, axis)
        G2 = nx.read_adjlist("tests/testcase1.adjlist")
        self.assertTrue(nx.is_isomorphic(G1, G2))

    def test_weighted_icosahedral_graph_generation(self):
        faceEdges = [
            ((1, 0), (1, 1)),  # c
            ((1, 0), (1, -1)),  # b
            ((1, 0), (0, 1)),  # b
            ((1, 0), (2, -1)),  # c
            ((1, 0), (0, 0)),  # a
            ((1, 0), (2, 0)),  # c
            ((1, 1), (2, 0)),  # b
            ((1, -1), (2, -1)),  # c
            ((2, 0), (2, 1)),  # a
            ((2, 0), (2, -1)),  # c
            ((2, 0), (3, -1)),  # c
            ((2, 0), (3, 0)),  # b
            ((2, -1), (2, -2)),  # b
            ((2, -1), (3, -2)),  # a
            ((2, -1), (3, -1)),  # b
        ]
        axis = ((0, 0), (2, 1), (3, -2))
        Ea = 1
        Eb = 2
        Ec = 3
        strength = [Ec, Eb, Eb, Ec, Ea, Ec, Eb, Ec, Ea, Ec, Ec, Eb, Eb, Ea, Eb]
        G1 = create_icosahedral_capsid_graph(faceEdges, axis, bond_strength=strength)
        G2 = nx.read_edgelist("tests/testcase2.edgelist")
        self.assertTrue(nx.is_isomorphic(G1, G2, edge_match=strength_edge_match))

    def test_create_icosahedral_graph_floating(self):
        P = icosahedral_patterns.PATTERN_6434
        h = 1
        k = 2
        [edges, Tx, Ty, Tscale] = P
        face_edges, axis = create_icosahedral_face_edges(edges, Tx, Ty, h, k)
        G1 = create_icosahedral_capsid_graph(face_edges, axis)
        self.assertEqual(20 * Tscale * (h * h + h * k + k * k), len(G1.nodes))
        self.assertTrue(
            nx.is_isomorphic(G1, nx.read_adjlist("tests/testcase3.adjlist"))
        )

    def test_create_cubic_graph(self):
        face = [
            ((0.25, -0.25), (-0.25, 0.25)),
            ((0.25, -0.25), (0.25, -0.75)),
            ((0.25, -0.25), (0.75, -0.25)),
            ((-0.25, 0.25), (-0.25, 0.75)),
            ((-0.25, 0.25), (-0.75, 0.25)),
        ]
        square_vertices = [
            (0.5, 0.5),
            (0.5, -0.5),
            (-0.5, -0.5),
            (-0.5, 0.5),
        ]
        w = [2, 1, 1, 1, 1]

        G1 = create_cubic_capsid_graph(face, square_vertices, w)
        G2 = nx.read_edgelist("tests/testcase4.edgelist")
        nx.is_isomorphic(G1, G2, edge_match=strength_edge_match)

    def test_texture_generator(self):
        from PIL import Image

        P = icosahedral_patterns.PATTERN_6434
        img = create_icosahedral_texture(P, 2, 1)
        img2 = Image.open("tests/test_texture.png")
        self.assertEqual(img.size, img2.size)
        self.assertEqual(img.mode, img2.mode)
        # test that pixels are the same
        img_matrix = np.asarray(img)
        img2_matrix = np.asarray(img2)
        self.assertTrue(np.all(img_matrix == img2_matrix))

    def test_cubic_generator(self):
        for P, testfile in [
            (cubic_patterns.AALS_24_PATTERN, "tests/AaLS_24.adjlist"),
            (cubic_patterns.AALS_48_PATTERN, "tests/AaLS_48.adjlist"),
            (cubic_patterns.AALS_60_PATTERN, "tests/AaLS_60.adjlist"),
        ]:
            [edges, Tx, Ty, face_side_edge] = P
            face_edges, face_square_vertices = create_cubic_face_edges(
                edges, Tx, Ty, face_side_edge
            )
            G = create_cubic_capsid_graph(face_edges, face_square_vertices)
            G2 = nx.read_adjlist(testfile)
            self.assertTrue(nx.is_isomorphic(G, G2))


if __name__ == "__main__":
    unittest.main()
