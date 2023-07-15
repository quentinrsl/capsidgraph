from typing import Dict, List, Tuple
import networkx as nx
from capsidgraph.util.types import Edge, Point, Vector2D
from capsidgraph.generator.face.util import is_under_line, extend


def rotate_point(P: Point, A: Point, cockwise: bool) -> Point:
    """
    Rotate a point P around a point A by 60° in the direction given by `cockwise`
    Coordinates are given in the "hexagonal" basis.

    Parameters
    ----------
    P : Point
        The point to rotate
    A : Point
        The center of the rotation
    cockwise : bool
        The direction of the rotation

    Returns
    -------
    Point
        The rotated point
    """
    if not cockwise:
        return (A[1] - P[1] + A[0], P[0] + P[1] - A[0])
    else:
        return (P[0] + P[1] - A[1], A[0] + A[1] - P[0])


def create_triangle(
    h: int, k: int, Tx: Vector2D, Ty: Vector2D
) -> Tuple[Point, Point, Point]:
    """
    Return the three points in the 2D plane corresponding to the vertices of the triangular face of the capsid.
    The coordinates are given in the "hexagonal" basis.

    Parameters
    ----------
    Tx : Vector2D
        The translation Vector2D between to hexagons in the tiling in the x direction
    Ty : Vector2D
        The translation Vector2D between to hexagons in the tiling in the y direction
    h : int
        The h Caspar Klug constant 
    k : int
        The k Caspar Klug constant

    Returns
    -------
    Tuple[Point, Point, Point]
        The three points of the triangular face in hexagnal coordinates.

    """
    A = (0, 0)
    B = (h * Tx[0] + k * Ty[0], h * Tx[1] + k * Ty[1])
    C = (B[0] + B[1], -B[0])  # Point B rotated around A
    return A, B, C


def is_in_triangle(P: Point, A: Point, B: Point, C: Point) -> int:
    """
    Return where is the point P relative to the triangle ABC.

    Parameters
    ----------
    P : Point
        The point to test
    A : Point
        The first vertex of the triangle
    B : Point
        The second vertex of the triangle
    C : Point
        The third vertex of the triangle

    Returns
    -------
    int
        -1 if P is outside the triangle, 0 if P is on the edge of the triangle, 1 if P is inside the triangle

    Notes
    -----
    All coordinates are given in the "hexagonal" basis.
    """
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


def extract_triangle(
    edges: List[Edge], h: int, k: int, Tx: Vector2D, Ty: Vector2D
) -> List[Edge]:
    """
    Extract all the edges that constitute of triangular face of the capsid.

    Parameters
    ----------
    edges : List[Edge]
        The edges of the tiling
    h : int
        The h Caspar Klug constant
    k : int
        The k Caspar Klug constant
    Tx : Vector2D
        The translation Vector2D between to hexagons in the tiling in the x direction
    Ty : Vector2D
        The translation Vector2D between to hexagons in the tiling in the y direction
    
    Returns
    -------
    List[Edge]
        The edges of the triangular face of the capsid

    Notes
    -----
    This returns all the segments of the tiling that have either : at least one extremity inside of the triangle or both extremities on the edge of the triangle.
    """
    res = []
    A, B, C = create_triangle(h, k, Tx, Ty)
    for P1, P2 in edges:
        if is_in_triangle(P1, A, B, C) + is_in_triangle(P2, A, B, C) >= 0:
            res.append((P1, P2))
    return res


def create_face_edges(
    edges: List[Edge], Tx: Vector2D, Ty: Vector2D, h: int, k: int
) -> Tuple[List[Edge], Tuple[Point, Point, Point]]:
    """
    Create the edges of the triangular face of a capsid, given a tiling pattern and the Caspar Klug constants.

    Parameters
    ----------
    edges : List[Edge]
        The edges of the one pattern of the tiling, in the "hexagonal" basis
    Tx : Vector2D
        The translation Vector2D between two patterns of the tiling in the x direction
    Ty : Vector2D
        The translation Vector2D between two patterns of the tiling in the y direction
    h : int
        The h Caspar Klug constant
    k : int
        The k Caspar Klug constant

    Returns
    -------
    Tuple[List[Edge], Tuple[Point, Point, Point]]
        The edges of the triangular face of the capsid, in the "hexagonal" basis, and the coordinates of the vertices of the triangular face
    """
    # Create a large enough tiling
    for i in range(h + k):
        edges = extend(edges, Tx)
    for i in range(h + k):
        edges = extend(edges, Ty)

    edges = extract_triangle(edges, h, k, Tx, Ty)
    return edges, create_triangle(h, k, Tx, Ty)
