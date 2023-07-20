from capsidgraph.util.types import Point, Vector2D, Edge
from typing import Tuple, List
from math import sqrt
from .util import extend, is_under_line

SQ3 = sqrt(3)


def rotate_point(P: Point, A: Point, cockwise: bool) -> Point:
    """
    Rotate a point P around a point A by 90° in the direction given by `cockwise`
    Coordinates are given in the canonical carthesian basis.

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
        return (A[1] - P[1] + A[0], -A[0] + P[0] + A[1])
    else:
        return (-A[1] + P[1] + A[0], A[0] - P[0] + A[1])


def create_square(edge: Vector2D) -> Tuple[Point, Point, Point, Point]:
    """
    Return the four points in the 2D plane corresponding to the vertices of the square face of the capsid.
    One of the points is at the origin.
    The rest of the points are deduced from the edge vector.

    Parameters
    ----------
    edge : Vector2D
        A vector representing one of the edges of the square face

    Returns
    -------
    Tuple[Point, Point, Point, Point]
        The four points of the sauare face in carthesian coordinates.

    """
    A = (0, 0)
    D = (edge[1], -edge[0])
    C = (edge[0] + edge[1], edge[1] - edge[0])
    return A, edge, C, D


def is_in_square(P: Point, A: Point, B: Point, C: Point, D: Point):
    """
    Return where is the point P relative to the square ABCD.

    Parameters
    ----------
    P : Point
        The point to test
    A : Point
        The first vertex of the square
    B : Point
        The second vertex of the square
    C : Point
        The third vertex of the square
    D : Point
        The fourth vertex of the square

    Returns
    -------
    int
        -1 if P is outside the square, 0 if P is on the edge of the square, 1 if P is inside the square
    """

    if (
        is_under_line(P, A, B) == -1
        or is_under_line(P, B, C) == -1
        or is_under_line(P, C, D) == -1
        or is_under_line(P, D, A) == -1
    ):
        return -1
    elif (
        is_under_line(P, A, B) == 0
        or is_under_line(P, B, C) == 0
        or is_under_line(P, C, D) == 0
        or is_under_line(P, D, A) == 0
    ):
        return 0
    else:
        return 1


def extract_square(edges: List[Point], square_edge: Vector2D):
    """
    Extract the square face from the list of edges.

    Parameters
    ----------
    edges : List[Point]
        The list of edges of the capsid
    square_edge : Vector2D
        The vector representing one of the edges of the square face
    
    Returns
    -------
    List[Point]
        The list of edges of the square face
    """
    res = []
    A, B, C, D = create_square(square_edge)
    for P1, P2 in edges:
        if is_in_square(P1, A, B, C, D) + is_in_square(P2, A, B, C, D) >= 0:
            res.append((P1, P2))
    return res


def create_face_edges(
    edges: List[Point], Tx: Vector2D, Ty: Vector2D, face_edge: Vector2D, repeat: int=1
) -> Tuple[List[Edge], Tuple[Point, Point, Point, Point]]:
    """
    Create the edges of the square face from a tiling pattern and a corner of the square face to use as the face of the cube.

    Parameters
    ----------
    edges : List[Point]
        The list of edges of the pattern
    Tx : Vector2D
        The first translation vector of the tiling pattern
    Ty : Vector2D
        The second translation vector of the tiling pattern
    face_edge : Vector2D
        The vector representing one of the edges of the square face
    repeat : int, optional
        The number of times to repeat the pattern, by default 10. Ensure that the pattern is large enough to cover the square face.

    Returns
    -------
    Tuple[List[Edge], Tuple[Point, Point, Point, Point]]
        The list of edges of the square face and the four vertices of the square face
    """

    for i in range(repeat):
        edges = extend(edges, Tx)
        edges = extend(edges, (-Tx[0], -Tx[1]))
    for i in range(repeat):
        edges = extend(edges, Ty)
        edges = extend(edges, (-Ty[0], -Ty[1]))

    edges = extract_square(edges, face_edge)
    return edges, create_square(face_edge)
