from multiprocessing import Pool
from matplotlib import pyplot as plt
import numpy as np

from capsidgraph.graph_generator import *
from capsidgraph.graph_analyser import *
from polyomavirus_generator import *

#Compute the fragmentation probability of the polyomavirus graph for different fractions of energy removed
#Save the results in results.txt and show the graph

#SIMULATION PARAMETERS
fb = 1 #fb = E_b / E_a
fc = 1 #fc = E_c / E_a

fragmentationType = "nodes"		#"nodes" or "edges", what to remove from the graph 
processCount = 4				#Number of process
iterationNumber = 1000000		#Nb d'iterations
maximumFragmentationEnergy=0.1	#Maximum value for the percolation energy
minimumFraglmentationEnergy=0.9	#Minmimum value for the percolation energy
pointNumber = 10				#Number of points to compute

a = 1 / (60 + 60 * fb + 90 * fc)
b = a * fb
c = a * fc

G = createPolyomavirusCapsid(a,b,c)
X = np.linspace(minimumFraglmentationEnergy,maximumFragmentationEnergy,pointNumber)


def worker(x):
	if(fragmentationType == "edges"):
		return getFragmentationWeightedProbabilityEdges(G,x, {"conditionType": 'fixedIterations', "iterations": iterationNumber},debug=True)[0]
	elif(fragmentationType == "nodes"):
		return getFragmentationWeightedProbabilityNodes(G,x, {"conditionType": 'fixedIterations', "iterations": iterationNumber},debug=True)[0]

if __name__ == '__main__':
	print(a,b,c)
	with Pool(processCount) as pool:
		Y=pool.map(worker,X)
	with open("results.txt",'a+') as f:
		f.write("fb : " + str(fb) + " fc : " + str(fc) + " n : " + str(iterationNumber) + "\n")
		f.write(str(list(X)) + " " +str(Y) + "\n")
	plt.plot(X,Y)
	plt.show()