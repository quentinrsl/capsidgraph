import unittest
from capsidgraph.generator import create_icosahedral_face_edges

class TestFaceGeneration(unittest.TestCase):
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

if __name__ == '__main__':
    unittest.main()