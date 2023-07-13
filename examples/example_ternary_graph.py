from multiprocessing import Pool
from capsidgraph.graph_generator import *
from capsidgraph.graph_analyser import *
import ternary
from polyomavirus_generator import createPolyomavirusCapsid

# Compute and display a ternary graph of the energy percolation threshold for different distribution of energy between different bond types.
# Settings
processes_number = 4
probaError = 0.05
nSteps = 8
size = 10
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
	#GRAPH
	## Boundary and Gridlines
	figure, tax = ternary.figure(scale=size)
	# Draw Boundary and Gridlines
	tax.boundary(linewidth=2.0)
	tax.get_axes().axis('off')
	tax.gridlines(color="black", multiple=5)
	tax.gridlines(color="blue", multiple=1, linewidth=0.5)

	# Set Axis labels and Title
	fontsize = 20
	# tax.set_title("Energy fragmentation threshold", fontsize=fontsize)
	tax.left_axis_label("$f_c$", fontsize=fontsize, offset=.15)
	tax.right_axis_label("$f_b$", fontsize=fontsize, offset=.15)
	tax.bottom_axis_label("$f_a$", fontsize=fontsize, offset=.15)

	# Set ticks
	# tax.ticks(axis='lbr', linewidth=1,multiple=5, offset=0.02)

	tax.heatmap(data, size,style="hexagon")

	#Draw lines
	p1 = (0,0,size)
	p2 = (2/7 * size, 2/7 * size, 3/7 * size)
	p3 = (2/5 * size, 0, 3/5 * size)
	tax.line(p1, p2, linewidth=2., color='red', linestyle="-")
	tax.line(p2, p3, linewidth=2., color='red', linestyle="-")
	tax.line(p1, p3, linewidth=2., color='red', linestyle="-")

	# Remove default Matplotlib Axes
	tax.clear_matplotlib_ticks()
	tax.show()

