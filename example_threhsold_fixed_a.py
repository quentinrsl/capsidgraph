from multiprocessing import Pool
import numpy as np
from dbcreds import *

from capsidlib.graph_generator import *
from capsidlib.graph_analyser import *
from polyomavirus_generator import *

#SIMULATION PARAMETERS
fragmentationType = "edges"

processCount = 1	#Number of processes
nSteps = 3			#Nb of iterations
pointNumber = 10	#Number of points to compute
probaError = 0.05
X = np.linspace(0.01,5/7 - 0.01,pointNumber)

def worker(x):
	a = (2/7)/60
	b = x/60
	c = (5/7-x)/90
	print(a,b,c)
	G = createPolyomavirusCapsid(a,b,c)
	if(fragmentationType == "edges"):
		th,N = getEnergyPercolationThresholdEdges(G,probaError,nSteps,debug=True,minIterations=1000)
	elif(fragmentationType == "nodes"):
		th,N = getEnergyPercolationThresholdNodes(G,probaError,nSteps,debug=True,minIterations=1000)
	return th

if __name__ == '__main__':
	with Pool(processCount) as pool:
		Y=pool.map(worker,X)
	print(X,Y)
	with open("results.txt",'a+') as f:
		f.write("#=====\n")
		f.write("#nSteps="+str(nSteps)+"perror="+str(probaError)+ " :\n")
		f.write("X="+str(X)+"\n")
		f.write("Y="+str(Y)+"\n")

