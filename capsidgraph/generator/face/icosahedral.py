from typing import Dict, List, Tuple
import networkx as nx
from capsidgraph.util.types import Edge, Point, Vector2D
from capsidgraph.generator.face.util import is_under_line, extend


def rotate_point(P: Point, A: Point, cockwise: bool) -> Point:
    if not cockwise:
        return (A[1] - P[1] + A[0], P[0] + P[1] - A[0])
    else:
        return (P[0] + P[1] - A[1], A[0] + A[1] - P[0])


# Return the three points in the 2D plane corresponding to the vertices of the triangular face of the capsid
# Tx et Ty are tuples of floats and represent the translation Vector2Ds between to hexagons in the tiling
# h and k are integer and represent the desired Caspar Klug configuration for the capsid
def create_triangle(
    h: int, k: int, Tx: Vector2D, Ty: Vector2D
) -> Tuple[Point, Point, Point]:
    A = (0, 0)
    B = (h * Tx[0] + k * Ty[0], h * Tx[1] + k * Ty[1])
    C = (B[0] + B[1], -B[0])  # Point B rotated around A
    return A, B, C


# Returns where is the point P relative to the triangle ABC
# A B C and P are tuples of floats
# -1 : outside the triangle
# 0 : sits on the edge of the triangle
# 1 : inside the triangle
def is_in_triangle(P: Point, A: Point, B: Point, C: Point) -> int:
    if (
        is_under_line(P, A, B) == -1
        or is_under_line(P, B, C) == -1
        or is_under_line(P, C, A) == -1
    ):
        return -1
    elif (
        is_under_line(P, A, B) == 0
        or is_under_line(P, B, C) == 0
        or is_under_line(P, C, A) == 0
    ):
        return 0
    else:
        return 1


# Extract a triangular face out of the tiling
# Edges is a list of tuples of floats and represent the edges of the tiling
# h and k are ntegers and represent the desired Caspar and Klug configuration for the capsid
# Tx et Ty are tuples of floats and represent the translation Vector2Ds between to hexagons in the tiling
# Returns all segments of the tiling that have either : at least one extremity inside of the triangle or both extremities on the edge of the triangle
def extract_triangle(
    edges: List[Edge], h: int, k: int, Tx: Vector2D, Ty: Vector2D
) -> List[Edge]:
    res = []
    A, B, C = create_triangle(h, k, Tx, Ty)
    for P1, P2 in edges:
        if is_in_triangle(P1, A, B, C) + is_in_triangle(P2, A, B, C) >= 0:
            res.append((P1, P2))
    return res


# Uses al the previously defined fucntions to create a list of edges corresponding to the triangular face of a capisd which tiling is made of a pattern given by the list of tuple edges,
# which needs to bez translated by Tx and Ty to create a tiling,
# and which Caspar Klug h and k numbers are the integers h and k
# Returns a list of tuples of floats and the coordinate of the vertices of the triangular face which those edges are part of
def create_face_edges(
    edges: List[Edge], Tx: Vector2D, Ty: Vector2D, h: int, k: int
) -> Tuple[List[Edge], Tuple[Point, Point, Point]]:
    # Create a large enough tiling
    for i in range(h + k):
        edges = extend(edges, Tx)
    for i in range(h + k):
        edges = extend(edges, Ty)

    edges = extract_triangle(edges, h, k, Tx, Ty)
    return edges, create_triangle(h, k, Tx, Ty)
