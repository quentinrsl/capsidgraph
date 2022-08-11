import random
import time
from typing import Dict, List, Tuple
import networkx as nx

#===========================
#non weighted graph analysis
#===========================

#returns the probability of the capsid beign fragmented if every nodes / edges are removed with probability p
#stopCondition is a dictionary that looks like this
# {
# 'conditionType': bisection / fixedIterations, 
# 'errorProbability': float between 0 and 1, only required for bisection, upper bound of the probability that the function will return a wrong result for the next bisection step (ie a value higer than 0.5 when the actual value is below 0.5 or the opposite, because in this stopCondition, we only care if the probability is above or below 0.5, not the actual vlaue).
# 'minIterations': minimum amount of iterations, only required for bisection
# 'maxIterations': maximum amount of iterations before giving up, only required for bisection, (to prevent unnecessary long computation, as the amount of ityeration to get above the lower bound of error is theoretically unbounded and gets bigger as we compute values closer to 0.5)
# 'iterations': number of iterations, only required for fixedIterations
# }
#fragType is "edges" or "nodes" and determines if edges or nodes will be removed
def getFragmentationProbability(G:nx.Graph,p:float,stopCondition:Dict,fragType:str,debug:bool=False)->int:
	start = time.time()
	fragmentedCount = 0
	pfrag = 0
	n = 0
	
	#Get simulation parameters
	conditionType = stopCondition["conditionType"]
	if(conditionType not in ["bisection", "fixedIterations"]):
		raise Exception("Invalid stop condition type provided")

	if(conditionType == "bisection"):
		errorProbability = stopCondition["errorProbability"]
		minIterations = stopCondition["minIterations"]
		maxIterations = stopCondition["maxIterations"]
		INV_PROBA=1/errorProbability
	#Condtion to satisfy to get under the upper bound of the probability of having a wrong result
	#Given by the Bienaymé-Tchebitchev inequality
	while (conditionType == "bisection" and ((4*n*((pfrag-0.5)**2) < INV_PROBA) or n<minIterations)) or (conditionType == "fixedIterations" and n < stopCondition['iterations']):
		#Give up if too many iterations
		if(conditionType == "bisection" and n > maxIterations):
			return (pfrag, True)
		if(n%1000 == 0 and debug):
			print("Progress : n=",n,"progress=", 100 * 4*n*((pfrag-0.5)**2) / INV_PROBA if conditionType == "bisection" else "")
		G_ = G.copy()
		if(fragType=="nodes"):
			#Remove each nodes with probability p
			for node in G.nodes:
				if(random.random() < p):
					G_.remove_node(node)
		elif(fragType=="edges"):
			#Remove each edge with probability p
			for a,b in G.edges:
				if(random.random() < p):
					G_.remove_edge(a,b)
		if (len(G_.nodes()) >0 and not nx.is_connected(G_)):
			fragmentedCount += 1
		n += 1
		#Compute probability
		pfrag = fragmentedCount / n
	if(debug):
		print("p=",p,"with n=",n,"got p(frag)=",pfrag,1000*(time.time()-start)/n,"ms/sim")
	return pfrag, False

#Returns the percolation threshold of a given graph with the probability of the result being wrong below errorProbability, as well as the number of bisection steps occured
#nSteps is an integer representing the number of bisection step
def getPercolationThreshold(G:nx.Graph,errorProbability:float,nSteps:int,fragType:str,minIterations:int,maxIterations:int)->Tuple[float,int]:
	eps = 1 - (1 - errorProbability) ** (1/nSteps)   #Compute the upper bond of error for one bisection step from the upper bond of making a mistake in the entire process
	a = 0
	b = 1
	n = 0
	while n < nSteps:
		m= (a+b)/2
		n += 1
		pfrag, aborted = getFragmentationProbability(G,m,{"conditionType":"bisection","errorProbability":eps,"minIterations":minIterations,"maxIterations":maxIterations},fragType)
		if aborted: 
			return m,n
		elif pfrag>0.5:
			b = m
		else:
			a = m
	return m,n

#Give an array A where A[i] is the number of fragments with size i encountered
#G is the graph to percolate on
#p is the probability of removal of each nodes/edges
#n is the number of iterations
#fragtype is either "nodes" or "edges" and determine if edges or nodes are to be removed 
def fragmentSizeDistribution(G:nx.Graph,p:float,n:int,fragType="nodes")->List[float]:
	fragmentSizes = {}
	maxSize = 0
	for i in range(n):
		G_ = G.copy()
		if(fragType=="nodes"):
		# On supprime chaque noeud avec une probability de p
			for node in G.nodes:
				if(random.random() < p):
					G_.remove_node(node)
		elif(fragType=="edges"):
			#On supprime chaque arrete avec une probability de p
			for a,b in G.edges:
				if(random.random() < p):
					G_.remove_edge(a,b)
		for component in nx.connected_components(G_):
			size = len(component)
			maxSize = max(size, maxSize)
			if not size in fragmentSizes:
				fragmentSizes[size] = 1
			else:
				fragmentSizes[size] += 1
	return [fragmentSizes[i] if i in fragmentSizes else 0 for i in range(maxSize - 1) ]

#===========================
# weighted graph analysis
#===========================

#returns the probability of the capsid beign fragmented if "perturbationEnergy" is removed from the capsid
#fragmentationEnergy: a float between 0 and 1 with 1 meaning all the capsid will be removed and 0 no energy will be removed
#stopCondition :
# {
# 'conditionType': bisection / fixedIterations, 
# 'errorProbability': float between 0 and 1, only required for bisection, upper bound of the probability that the bisection will give a wrong result for the next bisection step (ie a value higer than 0.5 when the actual value is below 0.5 or the opposite, because in this stopCondition, we only care if the probability is above or below 0.5, not the actual vlaue).
# 'minIterations': minimum amount of iterations, only required for bisection
# 'maxIterations': maximum amount of iterations before giving up, only required for bisection, (to prevent unnecessary long computation, as the amount of ityeration to get above the lower bound of error is theoretically unbounded and gets bigger as we compute values closer to 0.5)
# 'iterations': number of iterations, only required for fixedIterations
# }
#Returns a tuple (float, boolean), first value is the fragmentation energy, second value is whether or not the function stopped because it had reached the maximum amount of iterations
def getFragmentationWeightedProbabilityEdges(G:nx.Graph,fragmentationEnergy:float,stopCondition:Dict,debug:bool=False)->Tuple[float,bool]:
	#Get the attributes of the edges of the graph
	bondEnergy = nx.get_edge_attributes(G,'energy')
	edges = list(G.edges)
	#Compute edge probability weights
	weights = []
	for e in edges:
		weights.append(1/bondEnergy[e])
	minBondEnergy=1/max(weights)  #Energy of the weakest bond
	if(debug):
		print("STARTED E=",fragmentationEnergy,stopCondition)
	start = time.time()
	totalFragmented = 0 #Number of fragmented capids
	pfrag = 0   #Probability of fragmentation
	n = 0   #Amount of Monte Carlo simulations done

	#Get simulation parameters
	conditionType = stopCondition["conditionType"]
	if(conditionType not in ["bisection", "fixedIterations"]):
		raise Exception("Invalid stop condition type provided")
		
	if(conditionType == "bisection"):
		errorProbability = stopCondition["errorProbability"]
		minIterations = stopCondition["minIterations"]
		maxIterations = stopCondition["maxIterations"]
		INV_PROBA=1/errorProbability
	
	#Simulation
	#Check if stop condition is met
	while (conditionType == "bisection" and ((4*n*((pfrag-0.5)**2) < INV_PROBA) or n<minIterations)) or (conditionType == "fixedIterations" and n < stopCondition['iterations']):
		#Give up if too many iterations
		if(conditionType == "bisection" and n > maxIterations):
			return (pfrag, True)
		#Print debug statement
		if(n%1000 == 0 and debug):
			print("Progress : n=",n,"progress=", 100 * 4*n*((pfrag-0.5)**2) / INV_PROBA if conditionType == "bisection" else "")
		G_ = G.copy()
		remainingEdges = list(G_.edges)
		removedEdges = []
		remainigEdgesWeights = weights.copy()
		energy = fragmentationEnergy

		while energy > minBondEnergy :
			i = None
			#randomly pick a bond if there are any left, if the energy of the bond picked is higher than the percolationEnergy left, pick another one.
			#This loop has to stop because the energy is bigger than the minimum bond energy
			while (i == None or energy <= bondEnergy[remainingEdges[i]]) and len(remainingEdges) > 0:
				[i] = random.choices(list(range(len(remainingEdges))),weights=remainigEdgesWeights,k=1)

			#If a bond has been picked, remove it and make subsequent updates
			if i != None: 
				energy -= bondEnergy[remainingEdges[i]]
				removedEdges.append(remainingEdges[i])
				del remainingEdges[i]
				del remainigEdgesWeights[i]
				#update weakest bond energy
				if(len(remainigEdgesWeights) > 0):
					minBondEnergy=1/max(remainigEdgesWeights)
			else:
				#All edges have been removed
				break
		#Remove the edges from the graoh
		for (a,b) in removedEdges:
			G_.remove_edge(a,b)

		n += 1
		#Check if the graph is fragmented or not
		if (len(remainingEdges) > 0 and not nx.is_connected(G_)):
			totalFragmented += 1
		pfrag = totalFragmented/n
	#Another debug print
	if(debug):
		print("ENDED p(E=,",fragmentationEnergy,")=",pfrag,"time:",1000*(time.time()-start)/n,"ms/sim","N=",n)
	return (pfrag, False)

#Returns the percolation energy for which the probability of fragmentation is 0.5 of graph G, wit probability of beign wrong less than probaError and with a maximum amount of bisection step of nSteps
#The minIterations is the minimum amount of iterations before the next bisection step
#The maxIterations value is the maximum amount of iterations before the bisection process aborts
def getEnergyPercolationThresholdEdges(G:nx.Graph,probaError:float,nSteps:int,minIterations:int= 5000, maxIterations:int= 10000000,debug:bool=False)->Tuple[float,int]:
	eps = 1 - (1 - probaError) ** (1/nSteps)   #Compute the upper bond of error for one bisection step from the upper bond of making a mistake in the entire process
	a = 0
	b = 1
	n = 0
	while n < nSteps:
		m= (a+b)/2
		n += 1
		pfrag, aborted = getFragmentationWeightedProbabilityEdges(G,m,{"conditionType":"bisection","errorProbability":eps,"minIterations":minIterations,"maxIterations":maxIterations},debug=debug)
		if aborted: 
			return m,n
		elif pfrag > 0.5:
			b = m
		else:
			a = m
	return m,n

#intialize weights on nodes
#Create weight on nodes based on the sum of the energy of the bonds attached to it
#Ignore nodes that are in the blacklist (act if they were removed)hist
def getNodeAttributes(G:nx.Graph,node:int)->None:
		energy = 0
		for edge in G.edges(node):
			energy += G[edge[0]][edge[1]]['energy']
		return energy

#Same than with the edges but by percolating on nodes
#The energy of each node is the sum of the energy of the edges connected to it (the energy it would take to remove that node)
#We need to update the energy of the neigbors of a node that have been removed
def getFragmentationWeightedProbabilityNodes(G:nx.Graph,fragmentationEnergy:float,stopCondition:dict,debug:bool=False)->Tuple[float,bool]:
	#Initializing nodes attributes
	for node in G.nodes():
	# for node in G_.nodes():
		nodeEnergy = getNodeAttributes(G,node)
		G.nodes[node]["energy"] = nodeEnergy
	#Compute node weights
	#Compute edge probability weights
	weights = []
	for e in G.nodes:
		weights.append(1/G.nodes[e]["energy"])
	minNodeEnergy=1/max(weights)  #Energy of the weakest bond
	if(debug):
		print("STARTED NODES E=",fragmentationEnergy,stopCondition)

	start = time.time()
	totalFragmented = 0 #Number of fragmented capids
	pfrag = 0   #Probability of fragmentation
	n = 0   #Amount of Monte Carlo simulations

	#Get simulation parameters
	conditionType = stopCondition["conditionType"]
	if(conditionType not in ["bisection", "fixedIterations"]):
		raise Exception("Invalid stop condition type provided")
	if(conditionType == "bisection"):
		errorProbability = stopCondition["errorProbability"]
		minIterations = stopCondition["minIterations"]
		maxIterations = stopCondition["maxIterations"]
		INV_PROBA=1/errorProbability
	#Simulations
	while (conditionType == "bisection" and ((4*n*((pfrag-0.5)**2) < INV_PROBA) or n<minIterations)) or (conditionType == "fixedIterations" and n < stopCondition['iterations']):
		#Give up if too many iterations
		if(conditionType == "bisection" and n > maxIterations):
			return (pfrag, True)
		if(n%1000 == 0 and debug):
			print("Progress : n=",n,"progress=", 100 * 4*n*((pfrag-0.5)**2) / INV_PROBA if conditionType == "bisection" else "")
		G_ = G.copy()
		remainingNodes = list(G_.nodes)
		removedNodes = []
		remainingNodesWeights = weights.copy()
		energy = fragmentationEnergy

		minNodeEnergy=1/max(remainingNodesWeights)
		
		while energy + 1e-15 > minNodeEnergy:
			i = None
			while (i == None or energy + 1e-15 < G_.nodes[remainingNodes[i]]["energy"] ) and len(remainingNodes) > 0:
				[i] = random.choices(list(range(len(remainingNodes))),weights=remainingNodesWeights,k=1)
			if i != None:
				energy -= G_.nodes[remainingNodes[i]]["energy"]
				removedNodes.append(remainingNodes[i])
				# update the energy / probability weights of the neighbouring nodes
				removedNeighbours = []
				for (n1,n2,edgeAttributes) in G_.edges(remainingNodes[i],True):
					neighbour = n1 if n1 != remainingNodes[i] else n2
					if neighbour not in removedNodes:
						neighbourIndex = remainingNodes.index(neighbour)
						G_.nodes[neighbour]["energy"] -= edgeAttributes["energy"]

						#Put the neighbours to be removed in a separate list to keep the current node to remove at the index i in the remainingNode list
						if(G_.nodes[neighbour]["energy"] >1e-15):
							remainingNodesWeights[neighbourIndex] = abs(1/G_.nodes[neighbour]["energy"])
						else:
							removedNeighbours.append(remainingNodes[neighbourIndex])

				#remove node
				del remainingNodes[i]
				del remainingNodesWeights[i]
				#Remove 0 energy neighbours
				for neighbour in removedNeighbours:
					neighbourIndex = remainingNodes.index(neighbour)
					del remainingNodes[neighbourIndex]
					del remainingNodesWeights[neighbourIndex]
					removedNodes.append(neighbour)
				#update weakest bond energy
				if(len(remainingNodesWeights) > 0):
					minNodeEnergy=1/max(remainingNodesWeights)
			else:
				#All edges have been removed
				break

		for node in removedNodes:
			G_.remove_node(node)
		n += 1
		if (len(remainingNodes) > 0 and not nx.is_connected(G_)):
			totalFragmented += 1
		pfrag = totalFragmented/n
	if(debug):
		print("ENDED p(E=,",fragmentationEnergy,")=",pfrag,"time:",1000*(time.time()-start)/n,"ms/sim","N=",n)
	return (pfrag, False)    

#Same than with the edges but with nodes
def getEnergyPercolationThresholdNodes(G:nx.Graph,probaError:float,nSteps:int,minIterations:int= 5000, maxIterations:int= 10000000,debug:bool=False)->Tuple[float,int]:
	eps = 1 - (1 - probaError) ** (1/nSteps)   #Compute the upper bond of error for one bisection step from the upper bond of making a mistake in the entire process
	a = 0
	b = 1
	n = 0
	while n < nSteps:
		m= (a+b)/2
		n += 1
		pfrag, aborted = getFragmentationWeightedProbabilityNodes(G,m,{"conditionType":"bisection","errorProbability":eps,"minIterations":minIterations,"maxIterations":maxIterations},debug=debug)
		if aborted: 
			return m,n
		elif pfrag > 0.5:
			b = m
		else:
			a = m
	return m,n
