import unittest
import networkx as nx
import numpy as np

from capsidgraph.generator import create_icosahedral_face_edges
from capsidgraph.generator import create_icosahedral_capsid_graph
from capsidgraph.generator import create_cubic_capsid_graph
from capsidgraph.generator.texture.icosahedral import create_image
from math import sqrt

def energy_edge_match(e1,e2):
    return e1['energy']==e2['energy']
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
        face_edges, axis = create_icosahedral_face_edges(edges,Tx,Ty,h,k)
        self.assertEqual(face_edges,[((1, 0), (1, 1)), ((1, 1), (1, 2)), ((1, -1), (1, 0)), ((2, 0), (2, 1)), ((2, -1), (2, 0)), ((1, 0), (0, 1)), ((1, 1), (0, 2)), ((2, 0), (1, 1)), ((2, -1), (1, 0)), ((3, -1), (2, 0)), ((1, 0), (0, 0)), ((1, 1), (0, 1)), ((2, 0), (1, 0)), ((2, 1), (1, 1)), ((3, 0), (2, 0)), ((1, 0), (1, -1)), ((1, 1), (1, 0)), ((1, 2), (1, 1)), ((2, 0), (2, -1)), ((2, 1), (2, 0)), ((0, 1), (1, 0)), ((0, 2), (1, 1)), ((1, 0), (2, -1)), ((1, 1), (2, 0)), ((2, 0), (3, -1)), ((0, 0), (1, 0)), ((0, 1), (1, 1)), ((1, 0), (2, 0)), ((1, 1), (2, 1)), ((2, 0), (3, 0))])
        self.assertEqual(axis, ((0, 0), (1, 2), (3, -1)))
    
    def test_icosahedral_graph_generation(self):
        face_edges = [((1, 0), (1, 1)), ((1, 1), (1, 2)), ((1, -1), (1, 0)), ((2, 0), (2, 1)), ((2, -1), (2, 0)), ((1, 0), (0, 1)), ((1, 1), (0, 2)), ((2, 0), (1, 1)), ((2, -1), (1, 0)), ((3, -1), (2, 0)), ((1, 0), (0, 0)), ((1, 1), (0, 1)), ((2, 0), (1, 0)), ((2, 1), (1, 1)), ((3, 0), (2, 0)), ((1, 0), (1, -1)), ((1, 1), (1, 0)), ((1, 2), (1, 1)), ((2, 0), (2, -1)), ((2, 1), (2, 0)), ((0, 1), (1, 0)), ((0, 2), (1, 1)), ((1, 0), (2, -1)), ((1, 1), (2, 0)), ((2, 0), (3, -1)), ((0, 0), (1, 0)), ((0, 1), (1, 1)), ((1, 0), (2, 0)), ((1, 1), (2, 1)), ((2, 0), (3, 0))]
        axis =  ((0, 0), (1, 2), (3, -1))
        G1 = create_icosahedral_capsid_graph(face_edges,axis)
        G2 = nx.read_adjlist("tests/testcase1.adjlist")
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
        G2=nx.read_edgelist("tests/testcase2.edgelist")
        self.assertTrue(nx.is_isomorphic(G1,G2,edge_match=energy_edge_match))
    
    def test_create_icosahedral_graph_floating(self):
        IS3 = 1/sqrt(3)
        P = [
                [
                    ((1,0),(0,1)),
                    ((0,1),(-1,1)),
                    ((-1,1),(-1,0)),
                    ((-1,0),(0,-1)),
                    ((0,-1),(1,-1)),
                    ((1,-1),(1,0)),  
                    ((0,1),(-IS3,1+2*IS3)),
                    ((1,0),(1+IS3,IS3)),
                    ((0,1),(IS3,1+IS3)),
                    ((-1,1),(-1-IS3,1+2*IS3)),
                    ((1,0),(2*IS3+1,-IS3)),
                    ((1,-1),(2*IS3+1,-1-IS3)),
                ],
                (1+IS3,1+IS3),
                (-1-IS3, 2+2*IS3),
                3
        ]
        h = 1
        k = 2
        [edges,Tx,Ty,Tscale] = P
        face_edges, axis = create_icosahedral_face_edges(edges,Tx,Ty,h,k)
        G1 = create_icosahedral_capsid_graph(face_edges,axis)
        self.assertEqual(20 * Tscale * (h * h + h * k + k * k), len(G1.nodes))
        self.assertTrue(nx.is_isomorphic(G1,nx.read_adjlist("tests/testcase3.adjlist")))
    
    def test_create_cubic_graph(self):
        face = [
            ((0.25,-0.25),(-0.25,0.25)),
            ((0.25,-0.25),(0.25,-0.75)),
            ((0.25,-0.25),(0.75,-0.25)),
            ((-0.25,0.25),(-0.25,0.75)),
            ((-0.25,0.25),(-0.75,0.25))
        ]
        square_vertices = [
            (0.5,0.5),
            (0.5,-0.5),
            (-0.5,-0.5),
            (-0.5,0.5),
        ]
        w = [2,1,1,1,1]

        G1 = create_cubic_capsid_graph(face,square_vertices,w)
        G2 = nx.read_edgelist("tests/testcase4.edgelist")
        nx.is_isomorphic(G1,G2,edge_match=energy_edge_match)
    
    def test_texture_generator(self):
        IS3 = 1/sqrt(3)
        P = [
                [
                    ((1,0),(0,1)),
                    ((0,1),(-1,1)),
                    ((-1,1),(-1,0)),
                    ((-1,0),(0,-1)),
                    ((0,-1),(1,-1)),
                    ((1,-1),(1,0)),  
                    ((0,1),(-IS3,1+2*IS3)),
                    ((1,0),(1+IS3,IS3)),
                    ((0,1),(IS3,1+IS3)),
                    ((-1,1),(-1-IS3,1+2*IS3)),
                    ((1,0),(2*IS3+1,-IS3)),
                    ((1,-1),(2*IS3+1,-1-IS3)),
                ],
                (1+IS3,1+IS3),
                (-1-IS3, 2+2*IS3),
                3
        ]
        from PIL import Image
        img = create_image(P, 2, 1)
        img2 = Image.open("tests/test_texture.png")
        self.assertEqual(img.size,img2.size)
        self.assertEqual(img.mode,img2.mode)
        #test that pixels are the same
        img_matrix = np.asarray(img)
        img2_matrix = np.asarray(img2)
        self.assertTrue(np.all(img_matrix == img2_matrix))


if __name__ == '__main__':
    unittest.main()