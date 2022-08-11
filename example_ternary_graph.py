from multiprocessing import Pool
from capsidlib.graph_generator import *
from capsidlib.graph_analyser import *
import ternary
from polyomavirus_generator import createPolyomavirusCapsid
#SETTINGS
processes_number = 4
probaError = 0.05
nSteps = 8
size = 20
type = "nodes"

X = []
for i,j,k in ternary.helpers.simplex_iterator(size,False):
	X.append((i/size,j/size,k/size))

data = {}

def worker(param):
	(fa,fb,fc) = param
	a = fa/60
	b = fb/60
	c = fc/90
	if type == 'edges':
		th,N = getEnergyPercolationThresholdEdges(createPolyomavirusCapsid(a,b,c),probaError,nSteps,minIterations=1000)
	elif type == 'nodes':
		th,N = getEnergyPercolationThresholdNodes(createPolyomavirusCapsid(a,b,c),probaError,nSteps,minIterations=1000)
	data[(fa,fb,fc)] = th


if __name__ == '__main__':
	with Pool(processes_number) as pool:
		Y=pool.map(worker,X)
	with open("results.txt",'a+') as f:
		f.write("#=====\n")
		f.write("#nSteps="+str(nSteps)+"perror="+str(probaError)+ "size="+str(size)+ " :\n")
		f.write("data="+str(data)+"\n")

