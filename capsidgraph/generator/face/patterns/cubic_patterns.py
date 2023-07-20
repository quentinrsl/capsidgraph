from math import sqrt
from capsidgraph.generator.face.cubic import rotate_point

SQ3 = sqrt(3)
_AALS_24_POINTS = [
    (1 / 2, 1 / 2),
    (-1 / 2, 1 / 2),
    (-1 / 2, -1 / 2),
    (1 / 2, -1 / 2),
    (0, 1 / 2 + SQ3 / 2),
    (-SQ3 / 2, 1 + SQ3 / 2),
    (-SQ3 / 2 - 1, 1 + SQ3 / 2),
    (-SQ3 / 2 - 1 / 2, 1),
    (-SQ3 / 2 - 1 / 2, 0),
    (-SQ3 / 2 - 3 / 2, 0),
    (-SQ3 / 2 - 3 / 2, 1),
    (-SQ3 / 2 - 1, -SQ3 / 2),
    (-SQ3 / 2 - 1, -SQ3 / 2 - 1),
    (-1, -SQ3 / 2 - 1 / 2),
]

AALS_24_PATTERN = [
    [
        (_AALS_24_POINTS[0], _AALS_24_POINTS[1]),
        (_AALS_24_POINTS[1], _AALS_24_POINTS[2]),
        (_AALS_24_POINTS[2], _AALS_24_POINTS[3]),
        (_AALS_24_POINTS[3], _AALS_24_POINTS[0]),
        (_AALS_24_POINTS[0], _AALS_24_POINTS[4]),
        (_AALS_24_POINTS[4], _AALS_24_POINTS[5]),
        (_AALS_24_POINTS[5], _AALS_24_POINTS[6]),
        (_AALS_24_POINTS[6], _AALS_24_POINTS[7]),
        (_AALS_24_POINTS[7], _AALS_24_POINTS[1]),
        (_AALS_24_POINTS[7], _AALS_24_POINTS[8]),
        (_AALS_24_POINTS[8], _AALS_24_POINTS[9]),
        (_AALS_24_POINTS[9], _AALS_24_POINTS[10]),
        (_AALS_24_POINTS[10], _AALS_24_POINTS[7]),
        (_AALS_24_POINTS[8], _AALS_24_POINTS[1]),
        (_AALS_24_POINTS[11], _AALS_24_POINTS[8]),
        (_AALS_24_POINTS[12], _AALS_24_POINTS[11]),
        (_AALS_24_POINTS[13], _AALS_24_POINTS[12]),
        (_AALS_24_POINTS[2], _AALS_24_POINTS[13]),
    ],
    ((SQ3 + 1) / 2, -(SQ3 + 3) / 2),
    ((SQ3 + 3) / 2, (SQ3 + 1) / 2),
    (1 + SQ3 / 2, -1 / 2),
]

_AALS_48_POINTS = [
    (1 / 2, 1 / 2),
    (-1 / 2, 1 / 2),
    (0, 1 / 2 + SQ3 / 2),
    (1, 1 / 2 + SQ3 / 2),
    (0, 3 / 2 + SQ3 / 2),
    (1, 3 / 2 + SQ3 / 2),
    (1 / 2, 3 / 2 + SQ3),
    (1 / 2 + SQ3 / 2, 0),
    (1 + SQ3 / 2, SQ3 / 2),
    (1 + SQ3 / 2, SQ3 / 2 + 1),
]
_AALS_48_TILING = [
    (_AALS_48_POINTS[0], _AALS_48_POINTS[1]),
    (_AALS_48_POINTS[2], _AALS_48_POINTS[0]),
    (_AALS_48_POINTS[3], _AALS_48_POINTS[0]),
    (_AALS_48_POINTS[2], _AALS_48_POINTS[4]),
    (_AALS_48_POINTS[5], _AALS_48_POINTS[3]),
    (_AALS_48_POINTS[6], _AALS_48_POINTS[5]),
    (_AALS_48_POINTS[6], _AALS_48_POINTS[4]),
    (_AALS_48_POINTS[7], _AALS_48_POINTS[8]),
    (_AALS_48_POINTS[9], _AALS_48_POINTS[8]),
    (_AALS_48_POINTS[9], _AALS_48_POINTS[3]),
]

# Use rotational symetry of the tiling to complete it
for v, w in _AALS_48_TILING:
    # Rotated edge
    rv = rotate_point(v, (0, 0), False)
    rw = rotate_point(w, (0, 0), False)
    # Check if it already exists
    exists = False
    for v_, w_ in _AALS_48_TILING:
        if abs(rv[0] - v_[0]) < 1e-6 and abs(rv[1] - v_[1]) < 1e-6:
            if abs(rw[0] - w_[0]) < 1e-6 and abs(rw[1] - w_[1]) < 1e-6:
                exists = True
                break
        if abs(rv[0] - w_[0]) < 1e-6 and abs(rv[1] - w_[1]) < 1e-6:
            if abs(rw[0] - v_[0]) < 1e-6 and abs(rw[1] - v_[1]) < 1e-6:
                exists = True
                break
    if not exists:
        _AALS_48_TILING.append((rv, rw))

AALS_48_PATTERN = (
    _AALS_48_TILING,
    (2 + SQ3, -1),
    (1, 2 + SQ3),
    (3 / 2 + SQ3 / 2, 1 / 2 + SQ3 / 2),
)

AALS_60_PATTERN = (
    _AALS_48_TILING,
    (1.5 + SQ3, -2 - SQ3 / 2),
    (-2 - SQ3 / 2, -1.5 - SQ3),
    ( 7 / 4 + 3 * SQ3 / 4,(SQ3 - 1) / 4),
)
