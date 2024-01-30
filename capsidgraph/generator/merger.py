from typing import List, Dict, Tuple, Callable
import networkx as nx
from capsidgraph.util.types import Edge, Point

# Precision value to prevent floating point errors
eps = 1e-5

def _get_point_id(
    coordinates: Dict[int, List[Tuple[int, int, int]]],
    x: int,
    y: int,
    face_id: int | None = None,
) -> int:
    """
    Determine the id of a point given its coordinates in a given face if it exists, -1 otherwise.

    Parameters
    ----------
    coordinates : Dict[int, List[Tuple[int, int, int]]]
        The dictionnary of the points in the format `{id1:[(face1,x1,y1), (face2,x2,y2), ...], id2:...}`
    x : int
        The x coordinate of the point
    y : int
        The y coordinate of the point
    face_id : int
        The id of the face where its relative coordinates are (x,y) 
    """
    for id in coordinates:
        for comp_face_id, compx, compy in coordinates[id]:
            if (
                abs(compx - x) < eps
                and abs(compy - y) < eps
                and (face_id == None or comp_face_id == face_id)
            ):
                return id
    return -1


# Crée le graph d'une face a partir de segments,
# appartenants initalement à la face faceId
# Create a networkx graph of a single triangular face from a lattice pattern (edges)
# The ID of the nodes of the graph starts at startId
# Returns the graph and a dictionary saving to which node belong to which face (to each node id is linked an array of face Ids)
# This dictionary is initialized whith every nodes of the graph only belonging to the face faceId
def _create_face_graph(
    edges: List[Edge], face_id: int, start_id: int, bond_strength: List[float] = None
) -> Tuple[nx.Graph, Dict]:
    """
    Create the graph of a face given a list of segment.

    Parameters
    ----------
    Edges : List[Edges]
        Edges of a face to convert into a graph
    face_id : int
        The id of the face the graph represents
    bond_strength : List[float]
        List of the strength values for the edges, the i-th element of this list represent the strength of the i-th edge

    Returns
    -------
    Tuple[nx.Graph, Dict]
        The graph and the Dictionary containing the coordinates of each nodes
    """
    face_graph = nx.Graph()
    face_coordinates = {}
    # Create the nodes of the graph
    current_id = start_id
    for edge in edges:
        for x, y in edge:
            if _get_point_id(face_coordinates, x, y) == -1:
                face_graph.add_node(current_id)
                face_coordinates[current_id] = [(face_id, x, y)]
                current_id += 1

    # On crée the vertices of the graph
    for i in range(len(edges)):
        ((x1, y1), (x2, y2)) = edges[i]
        id1 = _get_point_id(face_coordinates, x1, y1)
        id2 = _get_point_id(face_coordinates, x2, y2)
        if bond_strength != None:
            w = bond_strength[i]
            face_graph.add_edge(id1, id2, strength=w)
        else:
            face_graph.add_edge(id1, id2)
    return face_graph, face_coordinates


def merge_faces(
    G: nx.Graph,
    coordinates: Dict[int, List[Tuple[int, int, int]]],
    rotate_point: Callable[[Point, Point, bool], Point],
    face_vertices: Tuple[Point, ...],
    id1: int,
    id2: int,
    vertex_id: int,
    clockwise: bool = True,
) -> nx.Graph:
    """
    Merge two faces in the graph by merging nodes belonging to both face id1 and id2.
    This is done by rotating face `id1` and rotating it around the verted `vertex_id` of the face polygon.
    It then matches the rotated points of `id1` with the points in face `id2` and merges points with corresponding coordinates.
    The `coordinates` Dict is updated in the process as well, by merging the array of matching points.
    Merged nodes now represent extremeties of edges belonging to multiple face at once.

    Parameters
    ----------
    G : nx.Graph
        The graph in which the faces will be merged
    coordinates : Dict[int, List[Tuple[int, int, int]]]
        The coordinate dictionnary matching node ids with the list of tuples in the format (face_id,x,y)
    rotate_point : Callable[[Point, Point, bool], Point]
        The function used to rotate points, it must take as parameters : the point to rotate, the center of rotation and the clockwise-direction of the rotation
    face_vertices : Tuple[Point, ...]
        A tuple of points representing the convex polygon of a "face" of the graph
    id1 : int
        id of the first face to merge
    id2 : int
        id of the second face to mergef
    vertex_id : int
        The index of the point of the face polygon around which rotate the face.
    clockwise : bool, optional
        Wether or not to rotate the face clockwise around the point, the value will be passed to the rotate_point function


    Returns
    -------
    nx.Graph
        The graph with the two faces merged
    """
    G_ = G
    nodes_to_remove = []
    # Get the points belonging to the face id1
    for id in coordinates:
        for face_id, x, y in coordinates[id]:
            if face_id == id1:
                rx, ry = rotate_point(
                    (x, y), face_vertices[vertex_id], clockwise
                )  # Rotate it by PI/3 in the right direction
                rid = _get_point_id(
                    coordinates, rx, ry, face_id=id2
                )  # Search for points in "id2" that would have the coordinate of this rotated point
                if rid != -1 and rid not in nodes_to_remove:  # We found a point
                    # Merge the points in the networkx graph as well as in the coordinate dict
                    coordinates[id] += coordinates[rid]
                    G_ = nx.contracted_nodes(G_, id, rid)
                    nodes_to_remove.append(rid)
    # Remove the dict entries after the search to not interfere with the loop
    for id in nodes_to_remove:
        # On supprime les entrées du dictionnaire dans un second temps pour ne pas perturber la boucle
        del coordinates[id]
    return G_
