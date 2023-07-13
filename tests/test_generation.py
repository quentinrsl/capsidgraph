import unittest
import networkx as nx

from capsidgraph.generator import create_icosahedral_face_edges
from capsidgraph.generator import create_icosahedral_capsid_graph

class TestGeneration(unittest.TestCase):
    def test_create_icosahedral_face_edges(self):
        pattern = [
            [
                ((0,0),(0,1)),
                ((0,0),(-1,1)),
                ((0,0),(-1,0)),
                ((0,0),(0,-1)),
                ((0,0),(1,-1)),
                ((0,0),(1,0)),
            ],
            (1, 0),
            (0,1),
            1
        ]
        h = 1
        k = 2
        [edges,Tx,Ty,Tscale] = pattern
        faceEdges, axis = create_icosahedral_face_edges(edges,Tx,Ty,h,k)
        self.assertEqual(faceEdges,[((1, 0), (1, 1)), ((1, 1), (1, 2)), ((1, -1), (1, 0)), ((2, 0), (2, 1)), ((2, -1), (2, 0)), ((1, 0), (0, 1)), ((1, 1), (0, 2)), ((2, 0), (1, 1)), ((2, -1), (1, 0)), ((3, -1), (2, 0)), ((1, 0), (0, 0)), ((1, 1), (0, 1)), ((2, 0), (1, 0)), ((2, 1), (1, 1)), ((3, 0), (2, 0)), ((1, 0), (1, -1)), ((1, 1), (1, 0)), ((1, 2), (1, 1)), ((2, 0), (2, -1)), ((2, 1), (2, 0)), ((0, 1), (1, 0)), ((0, 2), (1, 1)), ((1, 0), (2, -1)), ((1, 1), (2, 0)), ((2, 0), (3, -1)), ((0, 0), (1, 0)), ((0, 1), (1, 1)), ((1, 0), (2, 0)), ((1, 1), (2, 1)), ((2, 0), (3, 0))])
        self.assertEqual(axis, ((0, 0), (1, 2), (3, -1)))
    
    def test_icosahedral_graph_generation(self):
        face_edges = [((1, 0), (1, 1)), ((1, 1), (1, 2)), ((1, -1), (1, 0)), ((2, 0), (2, 1)), ((2, -1), (2, 0)), ((1, 0), (0, 1)), ((1, 1), (0, 2)), ((2, 0), (1, 1)), ((2, -1), (1, 0)), ((3, -1), (2, 0)), ((1, 0), (0, 0)), ((1, 1), (0, 1)), ((2, 0), (1, 0)), ((2, 1), (1, 1)), ((3, 0), (2, 0)), ((1, 0), (1, -1)), ((1, 1), (1, 0)), ((1, 2), (1, 1)), ((2, 0), (2, -1)), ((2, 1), (2, 0)), ((0, 1), (1, 0)), ((0, 2), (1, 1)), ((1, 0), (2, -1)), ((1, 1), (2, 0)), ((2, 0), (3, -1)), ((0, 0), (1, 0)), ((0, 1), (1, 1)), ((1, 0), (2, 0)), ((1, 1), (2, 1)), ((2, 0), (3, 0))]
        axis =  ((0, 0), (1, 2), (3, -1))
        G1 = create_icosahedral_capsid_graph(face_edges,axis)
        G2 = nx.read_adjlist("tests/testcase.adjlist")
        self.assertTrue(nx.is_isomorphic(G1,G2))
    
    def test_weighted_icosahedral_graph_generation(self):
        faceEdges = [((1,0),(1,1)), #c
            ((1,0),(1,-1)),#b
            ((1,0),(0,1)),#b
            ((1,0),(2,-1)),#c
            ((1,0),(0,0)),#a
            ((1,0),(2,0)),#c
            ((1,1),(2,0)),#b
            ((1,-1),(2,-1)),#c
            ((2,0),(2,1)),#a
            ((2,0),(2,-1)),#c
            ((2,0),(3,-1)),#c
            ((2,0),(3,0)),#b
            ((2,-1),(2,-2)),#b
            ((2,-1),(3,-2)),#a
            ((2,-1),(3,-1))#b
        ]
        axis = ((0, 0), (2, 1), (3, -2))
        Ea = 1
        Eb = 2
        Ec = 3
        energy =   [Ec,Eb,Eb,Ec,Ea,Ec,Eb,Ec,Ea,Ec,Ec,Eb,Eb,Ea,Eb]
        G1 = create_icosahedral_capsid_graph(faceEdges,axis,bond_strength=energy)
        G2=nx.read_edgelist("tests/testcase_weighted.edgelist")
        def edgematch(e1,e2):
            return e1['energy']==e2['energy']
        self.assertTrue(nx.is_isomorphic(G1,G2,edge_match=edgematch))

if __name__ == '__main__':
    unittest.main()