from capsidlib.graph_generator import *
from capsidlib.graph_analyser import *
from polyomavirus_generator import *
# Compute the energy percolation threshold for bond/edge removal, for given vales of fa,fb (The share of the total energy respectively give to types A and B bonds). 
# Settings :
fragmentationType = "edges"
nSteps = 8
probaError = 0.05
fa = 2/7
fb = 2/7
fc = 1-fa-fb

a = fa/60
b = fb/60
c = fc/90
print(a,b,c)
G = createPolyomavirusCapsid(a,b,c)
if(fragmentationType == "edges"):
    th,N = getEnergyPercolationThresholdEdges(G,probaError,nSteps,debug=True,minIterations=1000)
elif(fragmentationType == "nodes"):
    th,N = getEnergyPercolationThresholdNodes(G,probaError,nSteps,debug=True,minIterations=1000)
print("RESULT:",th,N)