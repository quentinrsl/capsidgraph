from capsidgraph.util.types import Point


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
