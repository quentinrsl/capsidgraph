import networkx as nx
import random

def probability(G: nx.Graph, settings: dict) -> nx.Graph:
	G_ = G.copy()
	p = settings.get("fragmentation_probability",0.5)
	fragmentation_type = settings.get("fragmentation_type","nodes")
	if(fragmentation_type=="nodes"):
		#Remove each nodes with probability p
		for node in G.nodes:
			if(random.random() < p):
				G_.remove_node(node)
	elif(fragmentation_type=="edges"):
		#Remove each edge with probability p
		for a,b in G.edges:
			if(random.random() < p):
				G_.remove_edge(a,b)
	return G_
