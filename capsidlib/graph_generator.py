from typing import Dict, List, Tuple
import networkx as nx

#Precision value to prevent floating point errors
eps = 1e-2

Point = Tuple[float,float]
Vector = Tuple[float,float]
Edge = Tuple[Point,Point]

#It is important to note that all coordonates in this algorithm are given in hexagonal coordinates to make rotations easier and the coordinates of most tilings easier to determine

#Duplicate all the edges and translate them by the vector T
#edges : list of tuples of floats, the list of edges
#T:  tuple of floats, the translation vector
#Returns the new list of edges
def extend(edges:List[Edge],T:Point)->List[Edge]:
	res = []
	for ((x1,y1),(x2,y2)) in edges:
		V1 = ((x1+T[0], y1+T[1]), (x2+T[0], y2+T[1]))
		V2 = ((x1-T[0], y1-T[1]), (x2-T[0], y2-T[1]))
		res.append(((x1,y1),(x2,y2)))
		if not V1 in edges:
			res.append(V1)
		if not V2 in edges:
			res.append(V2)
	return res

#Return the three points in the 2D plane corresponding to the vertices of the triangular face of the capsid 
#Tx et Ty are tuples of floats and represent the translation vectors between to hexagons in the tiling
#h and k are integer and represent the desired Caspar Klug configuration for the capsid
def createTriangle(h:int,k:int,Tx:Vector,Ty:Vector)->Tuple[Point,Point,Point]:
	A = (0,0)
	B = (h * Tx[0] + k * Ty[0], h * Tx[1] + k * Ty[1])
	C = (B[0] + B[1], -B[0])    #Point B rotated around A
	return A,B,C

#Returns whether or not the point P is under the line (AB)
#The order of the points A and B determine where "under" is
def isUnderLine(P:Point,A:Point,B:Point)->int:
	if (P[0]-A[0])*(B[1]-A[1]) > (P[1]-A[1])*(B[0]-A[0]): return 1
	if (P[0]-A[0])*(B[1]-A[1]) == (P[1]-A[1])*(B[0]-A[0]): return 0
	return -1

#Returns where is the point P relative to the triangle ABC
#A B C and P are tuples of floats
# -1 : outside the triangle
# 0 : sits on the edge of the triangle
# 1 : inside the triangle
def isInTriangle(P:Point,A:Point,B:Point,C:Point)->int:
	if isUnderLine(P,A,B) == -1 or isUnderLine(P,B,C) == -1 or isUnderLine(P,C,A) == -1:
		return -1
	elif isUnderLine(P,A,B) == 0 or isUnderLine(P,B,C) == 0 or isUnderLine(P,C,A) == 0:
		return 0
	else:
		return 1

#Extract a triangular face out of the tiling
#Edges is a list of tuples of floats and represent the edges of the tiling
#h and k are ntegers and represent the desired Caspar and Klug configuration for the capsid
#Tx et Ty are tuples of floats and represent the translation vectors between to hexagons in the tiling
#Returns all segments of the tiling that have either : at least one extremity inside of the triangle or both extremities on the edge of the triangle
def extractTriangle(edges:List[Edge],h:int,k:int,Tx:Vector,Ty:Vector)->List[Edge]:
	res = []
	A,B,C = createTriangle(h,k,Tx,Ty)
	for (P1,P2) in edges:
		if (isInTriangle(P1,A,B,C) + isInTriangle(P2,A,B,C) >= 0):
			res.append((P1,P2))
	return res


#Uses al the previously defined fucntions to create a list of edges corresponding to the triangular face of a capisd which tiling is made of a pattern given by the list of tuple edges, 
#which needs to bez translated by Tx and Ty to create a tiling, 
#and which Caspar Klug h and k numbers are the integers h and k
#Returns a list of tuples of floats and the coordinate of the vertices of the triangular face which those edges are part of
def createFaceEdges(edges:List[Edge], Tx:Vector, Ty:Vector, h:int, k:int)->Tuple[List[Edge],Tuple[Point,Point,Point]]:
	#Create a large enough tiling
	for i in range(h+k):
		edges = extend(edges,Tx)
	for i in range(h+k):
		edges = extend(edges,Ty)
	
	edges = extractTriangle(edges,h,k,Tx,Ty)
	return edges,createTriangle(h,k,Tx,Ty)
 
#Returns the coordinate of the point P rotated py PI/3 around the point A
def rotatePoint(P:Point,A:Point,cockwise:bool)->Point:
	if(not cockwise):
		return (
		A[1] - P[1] + A[0], 
		P[0] + P[1] - A[0]
		)
	else:
		return (
		P[0] + P[1] - A[1], 
		A[0] + A[1] - P[0]
		)

#Crée le graph d'une face a partir de segments, 
#appartenants initalement à la face faceId
#Create a networkx graph of a single triangular face from a lattice pattern (edges) 
#The ID of the nodes of the graph starts at startId 
#Returns the graph and a dictionary saving to which node belong to which face (to each node id is linked an array of face Ids) 
#This dictionary is initialized whith every nodes of the graph only belonging to the face faceId
def createFaceGraph(edges:List[Edge],faceId:int,startId:int,bondStrength:List[float]=None)->Tuple[nx.Graph,Dict]:
	faceGraph = nx.Graph()
	faceCoordonees = {}
	#determine if a point with given coordinates already exists, if so return its id, else return -1
	def getId(x,y):
		for id in faceCoordonees:
			for (_,compx,compy) in faceCoordonees[id]:
				if(abs(compx-x) < eps and abs(compy-y) < eps):
					return id
		return -1
	#Create the nodes of the graph
	id=startId
	for edge in edges:
		for (x,y) in edge:
			if getId(x,y) == -1:
				faceGraph.add_node(id)
				faceCoordonees[id] = [(faceId,x,y)]
				id += 1
	
	#On crée the vertices of the graph
	for i in range(len(edges)):
		((x1,y1),(x2,y2)) = edges[i]
		if(bondStrength != None):
			w=bondStrength[i]
		else:
			w=1
		id1 = getId(x1,y1)
		id2 = getId(x2,y2)
		faceGraph.add_edge(id1,id2,energy=w)
	return faceGraph,faceCoordonees


#Create the graph of the capsid
#faceEdges : an array of vertices corresponding to the edges in a traingular face of the capsid
#triangleVertices : tuple of three points corresponding to the three vertices of the triangular face overlayed on the lattice
def createCapsidGraph(faceEdges:List[Edge],triangleVertices:Tuple[Point,Point,Point],bondStrength:List[float]=None)->nx.Graph:
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
	def mergeFaces(id1:int,id2:int,vertexId:int,clockWise:bool=True)-> nx.Graph:
		G_ = G
		nodesToRemove = []
		#Get the points belonging to the face id1
		for id in coordinates:
			for (faceId,x,y) in coordinates[id]:
				if(faceId == id1):
					rx,ry = rotatePoint((x,y),triangleVertices[vertexId],clockWise) #Rotate it by PI/3 in the right direction
					rid = getId(id2,rx,ry)  #Search for points in "id2" that would have the coordinate of this rotated point
					if(rid != -1 and not rid in nodesToRemove): #We found a point
						#Merge the points in the networkx graph as well as in the coordinate dict
						coordinates[id] += coordinates[rid]
						G_ = nx.contracted_nodes(G_,id,rid)
						nodesToRemove.append(rid)
		#Remove the dict entries after the search to not interfere with the loop
		for id in nodesToRemove:
		#On supprime les entrées du dictionnaire dans un second temps pour ne pas perturber la boucle
			del coordinates[id]
		return G_

	#We then build the capsid by creating 20 triangular faces
	for faceId in range(20):
		face,faceCoo = createFaceGraph(faceEdges,faceId,nodeId,bondStrength=bondStrength)
		nodeId += len(face.nodes)
		G = nx.compose(G,face)
		coordinates.update(faceCoo)
	
	#Connect them according to the layout of an icosahedron
	#"Middle" row
	for i in range(5):
		G = mergeFaces(2*i,2*i+1,0,False)

	for i in range(5):
		G = mergeFaces(2*i+1,(2*i+2)%10,1,True)
	
	# "Top" and "Bottom" lines
	for i in range(10):
		G = mergeFaces(i, 10+i, (i+1)%2, i%2==0)
	
	#link faces together
	for i in range(10):
		#Replace the following line with the commeted ones for an incomplete build of the graph
		#Usefull for debuging
		# if(i%2==0):
		#     G = mergeFaces(10+i, 10 + ((i+2) % 10), 2, i%2==0)
		G = mergeFaces(10+i, 10 + ((i+2) % 10), 2, i%2==0)

	return G
