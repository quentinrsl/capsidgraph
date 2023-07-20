from capsidgraph.util.types import Edge, Point, Vector2D
from typing import Dict, List, Tuple
import networkx as nx


def extend(edges: List[Edge], T: Vector2D) -> List[Edge]:
    """
    Translate a list of edges by a vector T

    Parameters
    ----------
    edges : List[Edge]
        The list of edges to translate
    T : Vector2D
        The translation vector

    Returns
    -------
    List[Edge]
        The list of translated edges
    """
    res = []
    for (x1, y1), (x2, y2) in edges:
        V1 = ((x1 + T[0], y1 + T[1]), (x2 + T[0], y2 + T[1]))
        V2 = ((x1 - T[0], y1 - T[1]), (x2 - T[0], y2 - T[1]))
        res.append(((x1, y1), (x2, y2)))
        if V1 not in edges:
            res.append(V1)
        if V2 not in edges:
            res.append(V2)
    return res


def is_under_line(P: Point, A: Point, B: Point) -> int:
    """
    Determine whether a point P is under the line (AB) or not.
    The order of the points A and B determine where "under" is.

    Parameters
    ----------
    P : Point
        The point to locate
    A : Point
        The first point of the line
    B : Point
        The second point of the line

    Returns
    -------
    int
        1 if P is under (AB), 0 if P is on (AB), -1 if P is above (AB)
    """
    if (P[0] - A[0]) * (B[1] - A[1]) > (P[1] - A[1]) * (B[0] - A[0]):
        return 1
    if (P[0] - A[0]) * (B[1] - A[1]) == (P[1] - A[1]) * (B[0] - A[0]):
        return 0
    return -1
