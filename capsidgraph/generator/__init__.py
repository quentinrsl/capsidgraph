import networkx as nx
from capsidgraph.util.types import *
from typing import List, Tuple, Dict

from .face.icosahedral import rotate_point as _rotate_icosahedral_point
from .face.icosahedral import create_face_edges as create_icosahedral_face_edges
from .face.patterns import icosahedral_patterns

from .face.cubic import rotate_point as _rotate_cubic_point
from .face.cubic import create_face_edges as create_cubic_face_edges
from .face.patterns import cubic_patterns

from .merger import merge_faces, _create_face_graph

from .exceptions import GraphCreationException

from .texture.icosahedral import create_texture as create_icosahedral_texture


def create_icosahedral_capsid_graph(
    face_edges: List[Edge],
    triangle_vertices: Tuple[Point, Point, Point],
    bond_strength: List[float] | None = None,
) -> nx.Graph:
    if bond_strength != None and 0 in bond_strength:
        raise GraphCreationException("Bond cannot have 0 strength")
    G = nx.Graph()
    coordinates = {}
    node_id = 0
    # We then build the capsid by creating 20 triangular faces
    for face_id in range(20):
        face, face_coordinates = _create_face_graph(
            face_edges, face_id, node_id, bond_strength=bond_strength
        )
        node_id += len(face.nodes)
        G = nx.compose(G, face)
        coordinates.update(face_coordinates)

    # Connect them according to the layout of an icosahedron
    # "Middle" row
    for i in range(5):
        G = merge_faces(
            G,
            coordinates,
            _rotate_icosahedral_point,
            triangle_vertices,
            2 * i,
            2 * i + 1,
            0,
            False,
        )
    for i in range(5):
        G = merge_faces(
            G,
            coordinates,
            _rotate_icosahedral_point,
            triangle_vertices,
            2 * i + 1,
            (2 * i + 2) % 10,
            1,
            True,
        )
    # "Top" and "Bottom" lines
    for i in range(10):
        G = merge_faces(
            G,
            coordinates,
            _rotate_icosahedral_point,
            triangle_vertices,
            i,
            10 + i,
            (i + 1) % 2,
            i % 2 == 0,
        )
    # link faces together
    for i in range(10):
        # Replace the following line with the commeted ones for an incomplete build of the graph
        G = merge_faces(
            G,
            coordinates,
            _rotate_icosahedral_point,
            triangle_vertices,
            10 + i,
            10 + ((i + 2) % 10),
            2,
            i % 2 == 0,
        )
    return G


def create_cubic_capsid_graph(
    face_edges: List[Edge],
    square_vertices: Tuple[Point, Point, Point, Point],
    bond_strength: List[float] | None = None,
) -> nx.Graph:
    if bond_strength != None and 0 in bond_strength:
        raise GraphCreationException("Bond cannot have 0 strength")
    G = nx.Graph()
    coordinates = {}
    node_id = 0
    for face_id in range(6):
        face, face_coordinates = _create_face_graph(
            face_edges, face_id, node_id, bond_strength=bond_strength
        )
        node_id += len(face.nodes)
        G = nx.compose(G, face)
        coordinates.update(face_coordinates)

    G = merge_faces(G, coordinates, _rotate_cubic_point, square_vertices, 0, 1, 0, True)
    G = merge_faces(
        G, coordinates, _rotate_cubic_point, square_vertices, 1, 2, 2, False
    )
    G = merge_faces(G, coordinates, _rotate_cubic_point, square_vertices, 2, 3, 0, True)
    G = merge_faces(
        G, coordinates, _rotate_cubic_point, square_vertices, 3, 0, 2, False
    )
    G = merge_faces(G, coordinates, _rotate_cubic_point, square_vertices, 0, 4, 1, True)
    G = merge_faces(
        G, coordinates, _rotate_cubic_point, square_vertices, 4, 1, 2, False
    )
    G = merge_faces(G, coordinates, _rotate_cubic_point, square_vertices, 4, 2, 2, True)
    G = merge_faces(G, coordinates, _rotate_cubic_point, square_vertices, 4, 3, 3, True)
    G = merge_faces(G, coordinates, _rotate_cubic_point, square_vertices, 0, 5, 3, True)
    G = merge_faces(G, coordinates, _rotate_cubic_point, square_vertices, 5, 1, 1, True)
    G = merge_faces(G, coordinates, _rotate_cubic_point, square_vertices, 5, 2, 0, True)
    G = merge_faces(
        G, coordinates, _rotate_cubic_point, square_vertices, 5, 3, 0, False
    )

    return G
