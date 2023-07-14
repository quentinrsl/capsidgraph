from capsidgraph.util.types import Point


def rotate_point(P: Point, A: Point, cockwise: bool) -> Point:
    if not cockwise:
        return (A[1] - P[1] + A[0], -A[0] + P[0] + A[1])
    else:
        return (-A[1] + P[1] + A[0], A[0] - P[0] + A[1])
