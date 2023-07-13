from capsidgraph.util.types import *
from typing import Dict, List, Tuple
import networkx as nx


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

#Returns whether or not the point P is under the line (AB)
#The order of the points A and B determine where "under" is
def is_under_line(P:Point,A:Point,B:Point)->int:
	if (P[0]-A[0])*(B[1]-A[1]) > (P[1]-A[1])*(B[0]-A[0]): return 1
	if (P[0]-A[0])*(B[1]-A[1]) == (P[1]-A[1])*(B[0]-A[0]): return 0
	return -1