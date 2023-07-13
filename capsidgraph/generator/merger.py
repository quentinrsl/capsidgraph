from typing import List,Dict,Tuple,Callable
import networkx as nx
from util.types import *

#Precision value to prevent floating point errors
eps = 1e-5

#Returns the id of a given point if it exists, -1 otherwise
def get_point_id(coordinates:Dict[int,List[Tuple[int,int,int]]],x:int,y:int,face_id:int|None=None)->int:
    for id in coordinates:
        for (comp_face_id,compx,compy) in coordinates[id]:
            if(abs(compx-x) < eps and abs(compy-y) < eps and (comp_face_id == None or comp_face_id == face_id)):
                return id
    return -1
#Crée le graph d'une face a partir de segments, 
#appartenants initalement à la face faceId
#Create a networkx graph of a single triangular face from a lattice pattern (edges) 
#The ID of the nodes of the graph starts at startId 
#Returns the graph and a dictionary saving to which node belong to which face (to each node id is linked an array of face Ids) 
#This dictionary is initialized whith every nodes of the graph only belonging to the face faceId
def create_face_graph(edges:List[Edge],face_id:int,start_id:int,bond_strength:List[float]=None)->Tuple[nx.Graph,Dict]:
	face_graph = nx.Graph()
	face_coordinates = {}
	#Create the nodes of the graph
	current_id=start_id
	for edge in edges:
		for (x,y) in edge:
			if get_point_id(face_coordinates,x,y) == -1:
				face_graph.add_node(current_id)
				face_coordinates[current_id] = [(face_id,x,y)]
				current_id += 1
	
	#On crée the vertices of the graph
	for i in range(len(edges)):
		((x1,y1),(x2,y2)) = edges[i]
		if(bond_strength != None):
			w=bond_strength[i]
		else:
			w=1
		id1 = get_point_id(face_coordinates,x1,y1)
		id2 = get_point_id(face_coordinates,x2,y2)
		face_graph.add_edge(id1,id2,energy=w)
	return face_graph,face_coordinates



#Merge two faces in the graph
#This function takes the face with id "id1", rotated it clockwise / anticlocksie around the vertex of the triangle with vertexId (0,1 or 2)
#Then matches the points of this rotated face with the face with id "id2" and merges the point that have the same coordinates
#This update the coordinate dictionary as well by merging the "belonging arrays" of the two old nodes into the new node
#That way a node can represent the extremity of different edges belonging to different faces 
def merge_faces(G:nx.Graph,coordinates:Dict[int,List[Tuple[int,int,int]]],rotate_point:Callable[[Point,Point,bool],Point],face_vertices:Tuple[Point, ...],id1:int,id2:int,vertex_id:int,clockwise:bool=True)-> nx.Graph:
    G_ = G
    nodes_to_remove = []
    #Get the points belonging to the face id1
    for id in coordinates:
        for (face_id,x,y) in coordinates[id]:
            if(face_id == id1):
                rx,ry = rotate_point((x,y),face_vertices[vertex_id],clockwise) #Rotate it by PI/3 in the right direction
                rid = get_point_id(coordinates,rx,ry,face_id=id2)  #Search for points in "id2" that would have the coordinate of this rotated point
                if(rid != -1 and rid not in nodes_to_remove): #We found a point
                    #Merge the points in the networkx graph as well as in the coordinate dict
                    coordinates[id] += coordinates[rid]
                    G_ = nx.contracted_nodes(G_,id,rid)
                    nodes_to_remove.append(rid)
    #Remove the dict entries after the search to not interfere with the loop
    for id in nodes_to_remove:
    #On supprime les entrées du dictionnaire dans un second temps pour ne pas perturber la boucle
        del coordinates[id]
    return G_