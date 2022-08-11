from capsidlib.graph_generator import *
from capsidlib.graph_analyser import *
from polyomavirus_generator import *

#SIMULATION PARAMETERS
fragmentationType = "nodes"
nSteps = 8
probaError = 0.05
fa = 2/7
fb = 2/7
fc = 3/7

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