from capsidgraph.graph_generator import createFaceGraph
from typing import Dict, List, Tuple
import networkx as nx

#Precision value to prevent floating point errors
eps = 1e-2

Point = Tuple[float,float]
Vector = Tuple[float,float]
Edge = Tuple[Point,Point]

#Returns the coordinate of the point P rotated py 90° around the point A
def rotatePoint(P:Point,A:Point,cockwise:bool)->Point:
	if(not cockwise):
		return (
		A[1] - P[1] + A[0], 
		- A[0] + P[0] + A[1]
		)
	else:
		return (
		- A[1] + P[1] + A[0], 
		A[0] - P[0] + A[1]
		)

#Create the graph
#faceEdges : an array of vertices corresponding to the edges in a square face
#squareVertices: tuple of four points corresponding to the four vertices of the triangular face overlayed on the lattice
def createCapsidGraph(faceEdges:List[Edge],  squareVertices:Tuple[Point,Point,Point,Point], bondStrength:List[float]|None=None)->nx.Graph:
	if bondStrength != None and 0 in bondStrength:
		raise Exception("Bond cannot have 0 energy")
	
	G = nx.Graph()
	coordinates = {}
	nodeId = 0
	#Returns the id of a given point if it exists, -1 otherwise
	def getId(faceId:int,x:int,y:int)->int:
		for id in coordinates:
			for (compFaceId,compx,compy) in coordinates[id]:
				if(abs(compx-x) < eps and abs(compy-y) < eps and compFaceId == faceId):
					return id
		return -1

	#Merge two faces in the graph
	#This function takes the face with id "id1", rotated it clockwise / anticlocksie around the vertex of the triangle with vertexId (0,1 or 2)
	#Then matches the points of this rotated face with the face with id "id2" and merges the point that have the same coordinates
	#This update the coordinate dictionary as well by merging the "belonging arrays" of the two old nodes into the new node
	#That way a node can represent the extremity of different edges belonging to different faces 
	def mergeFaces(G:nx.graph,id1:int,id2:int,vertexId:int,clockWise:bool=True)-> nx.Graph:
		nodesToRemove = []
		#Get the points belonging to the face id1
		for id in coordinates:
			for (faceId,x,y) in coordinates[id]:
				if(faceId == id1):
					rx,ry = rotatePoint((x,y),squareVertices[vertexId],clockWise) #Rotate it by 90 degrees in the right direction
					rid = getId(id2,rx,ry)  #Search for points in "id2" that would have the coordinate of this rotated point
					if(rid != -1 and not rid in nodesToRemove): #We found a point
						#Merge the points in the networkx graph as well as in the coordinate dict
						coordinates[id] += coordinates[rid]
						G = nx.contracted_nodes(G,id,rid)
						nodesToRemove.append(rid)
		#Remove the dict entries after the search to not interfere with the loop
		for id in nodesToRemove:
		#On supprime les entrées du dictionnaire dans un second temps pour ne pas perturber la boucle
			del coordinates[id]
		return G

	#We then build the capsid by creating 20 triangular faces
	for faceId in range(6):
		face,faceCoo = createFaceGraph(faceEdges,faceId,nodeId,bondStrength=bondStrength)
		nodeId += len(face.nodes)
		G = nx.compose(G,face)
		coordinates.update(faceCoo)
	
	#Connect them according to the layout of an icosahedron

	G = mergeFaces(G,0,1,0,True)	
	G = mergeFaces(G,1,2,2,False)	
	G = mergeFaces(G,2,3,0,True)	
	G = mergeFaces(G,3,0,2,False)	

	G = mergeFaces(G,0,4,1,True)
	G = mergeFaces(G,4,1,2,False)
	G = mergeFaces(G,4,2,2,True)
	G = mergeFaces(G,4,3,3,True)
	
	G = mergeFaces(G,0,5,3,True)
	G = mergeFaces(G,5,1,1,True)
	G = mergeFaces(G,5,2,0,True)
	G = mergeFaces(G,5,3,0,False)

	return G