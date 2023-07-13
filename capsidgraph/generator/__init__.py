import networkx as nx
from capsidgraph.util.types import *
from typing import List, Tuple, Dict

from  .face.icosahedral import create_face_edges as create_icosahedral_face_edges
from .face.icosahedral import rotate_point as rotate_icosahedral_point
from .merger import merge_faces,create_face_graph

def create_icosahedral_capsid_graph(faceEdges:List[Edge],triangleVertices:Tuple[Point,Point,Point],bond_strength:List[float] | None=None)->nx.Graph:
	if bond_strength != None and 0 in bond_strength:
		raise Exception("Bond cannot have 0 energy")
	G = nx.Graph()
	coordinates = {}
	node_id = 0
	#We then build the capsid by creating 20 triangular faces
	for face_id in range(20):
		face,face_coordinates = create_face_graph(faceEdges,face_id,node_id,bond_strength=bond_strength)
		node_id += len(face.nodes)
		G = nx.compose(G,face)
		coordinates.update(face_coordinates)
	
	#Connect them according to the layout of an icosahedron
	#"Middle" row
	for i in range(5):
		G = merge_faces(G,coordinates, rotate_icosahedral_point,triangleVertices,2*i,2*i+1,0,False)

	for i in range(5):
		G = merge_faces(G,coordinates, rotate_icosahedral_point,triangleVertices,2*i+1,(2*i+2)%10,1,True)
	
	# "Top" and "Bottom" lines
	for i in range(10):
		G = merge_faces(G,coordinates, rotate_icosahedral_point,triangleVertices,i, 10+i, (i+1)%2, i%2==0)
	
	#link faces together
	for i in range(10):
		#Replace the following line with the commeted ones for an incomplete build of the graph
		#Usefull for debuging
		# if(i%2==0):
		#     G = mergeFaces(10+i, 10 + ((i+2) % 10), 2, i%2==0)
		G = merge_faces(G,coordinates, rotate_icosahedral_point,triangleVertices,10+i, 10 + ((i+2) % 10), 2, i%2==0)
	return G


def createCubicCapsidGraph(faceEdges:List[Edge],  squareVertices:Tuple[Point,Point,Point,Point], bondStrength:List[float]|None=None)->nx.Graph:
    pass