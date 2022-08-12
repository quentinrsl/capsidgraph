import matplotlib.pyplot as plt
from math import sqrt

from capsidlib.graph_generator import *
from capsidlib.graph_analyser import *

pattern666 = [
	[
		#Pattern to repeat
		((1,0),(0,1)),
		((0,1),(-1,1)),
		((-1,1),(-1,0)),
		((-1,0),(0,-1)),
		((0,-1),(1,-1)),
		((1,-1),(1,0))
	],
	#Translation vector to create the lattice
	(1, 1),
	(-1,2),
	#Scaling factor to compute triangulation number
	3
]

pattern666dual = [
	[
		((0,0),(0,1)),
		((0,0),(-1,1)),
		((0,0),(-1,0)),
		((0,0),(0,-1)),
		((0,0),(1,-1)),
		((0,0),(1,0)),
	],
	(1, 0),
	(0,1),
	1
]

pattern6363 = [
		[
			((1,0),(0,1)),
			((0,1),(-1,1)),
			((-1,1),(-1,0)),
			((-1,0),(0,-1)),
			((0,-1),(1,-1)),
			((1,-1),(1,0))
		],
		(2,0),
		(0,2),
		2
]

pattern3336 = [
		[
				((1,0),(0,1)),
				((0,1),(-1,1)),
				((-1,1),(-1,0)),
				((-1,0),(0,-1)),
				((0,-1),(1,-1)),
				((1,-1),(1,0)),

				((0,1),(0,2)),
				((0,2),(-1,2)),
				((0,1),(-1,2)),
				((-1,1),(-1,2)),
				((0,1),(1,1)),
				((0,2),(1,1)),
				((1,0),(1,1)),
				((1,0),(2,0)),
				((1,1),(2,0)),
				((2,-1),(2,0)),
				((2,-1),(1,0)),
				((1,-1),(2,-1)),
		],
		(2,1),
		(-1,3),
		3
]


IS3 = 1/sqrt(3)
pattern6434 = [
		[
				((1,0),(0,1)),
				((0,1),(-1,1)),
				((-1,1),(-1,0)),
				((-1,0),(0,-1)),
				((0,-1),(1,-1)),
				((1,-1),(1,0)),  
				((0,1),(-IS3,1+2*IS3)),
				((1,0),(1+IS3,IS3)),
				((0,1),(IS3,1+IS3)),
				((-1,1),(-1-IS3,1+2*IS3)),
				((1,0),(2*IS3+1,-IS3)),
				((1,-1),(2*IS3+1,-1-IS3)),
		],
		(1+IS3,1+IS3),
		(-1-IS3, 2+2*IS3),
		3
]

# Generation of a graph with T=7 and a triangular lattice
P = pattern666dual
h = 1
k = 2
[edges,Tx,Ty,Tscale] = P
faceEdges, axis = createFaceEdges(edges,Tx,Ty,h,k)
print("T=",Tscale * (h * h + h * k + k * k))
G = createCapsidGraph(faceEdges,axis)
nx.draw_kamada_kawai(G,node_size=5)
plt.show()