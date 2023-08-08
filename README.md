# Capsidgraph
Capsidgraph is a python library for generating and studying the graphs representing the protein subunit interaction networks of viral capsids.
# Requirements
* Python 3.10 or higher
* The `networkx` library is used to handle operations on graphs.
* The `matplotlib` and `numpy` libraries are used to create figures.
* To use the capsidgraph library, you need to add the capsidgraph folder to your PYTHONPATH environment variable. This can be done by adding the following line to your .bashrc file (replace the path with the path to the capsidgraph folder on your computer):
```bash
export PYTHONPATH=$PYTHONPATH:/path/to/capsidgraph
```

# Graph analysis
The `capsidgraph.analyser` module provides tools to compute information on the bond strength of different capsid graphs. This module requires networkx graphs with weighted nodes and can give insight into the ability of a capsid to resist both the removal of bonds (edges in the graph) and capsomers, i.e. protein units corresponding to the building blocks of the capsid (nodes in the graph).

## Unweighted graphs
### Fragmentation probability
For a graph that does not have weighted edges, the functions `get_fragmentation_probability_random_node_removal` and `get_fragmentation_probability_random_edge_removal` approximate the probability $p_f$ that the graph $G$ will fragment when a node or edge is removed with probability $f_r$. The value of $p_f$ is determined using a Monte Carlo method: for each simulation we determine whether or not to remove a node/edge with the built-in `random` function, and then check the connectivity of the graph using the networkx function `nx.is_connected`.

The function `probability_fragment` implements the random edge or node removal process. It takes as argument the graph to fragment and a Dict contaning two entries : 
- `fragmentation`: float between 0 and 1, probability of removal
- `fragmentation_type` : a string equals to "nodes" or "edges", determine wether to remove nodes or edges from the graph.


### Percolation threshold
The functions `get_fragmentation_probability_threshold_node` and `get_fragmentation_probability_threshold_edge` use a bisection method to approximate the threshold probability of removal for which the probability of fragmentation is 0.5, i.e the percolation threshold. 

The parameter `error_probability` provides an upper bound for the probability of the predicted value being correct. This routine uses `get_fragmentation_probability` to determine whether or not the probability of fragmentation for a given probability of removal is above or below 0.5 with a high enough probability.
To ensure that the probability of the value being incorrect is below `error_probability`, we stop a bisection step if the following inequality is true $4N|\frac{S_N}{N}-0.5|^2>\frac{1}{\epsilon}$ (where $\epsilon$ is the upper bond of the probability of error, $S_N$ the number of simulations that led to a fragmented graph and $N$ the number of simulations, for more details on this see the paper)

Minimum and maximum numbers of iterations for each bisection step can be provided as well. If the maximum number of simulations is reached before the probability condition is met, the bisection algorithm stops.

## Weighted graphs

### Removing edges
The energy of each bond is stored as the `energy` attribute of the corresponding edge in the graph. When removing edges in the graph we take into account this energy by attributing probability weights to each edge. This weight is defined as being proportional to the inverse of the energy. 
In order to remove a certain “energy amount” from the capsid, we define the notion of fragmentation probability by repeating following procedure:
* We first randomly pick an edge from the graph using the build it `choices` function from the `random` library.
* We then remove the energy of this bond from the total, thus computing the amount of energy that is left to remove.
* We remove the corresponding bond from the list of remaining bonds.
* We repeat this process as long as there is a bond in the graph that has less energy than the amount we have yet to remove.
* Once we are no longer able to remove edges, we determine whether or not the graph is fragmented using the `networkx` method for both removing edges and checking the connectivity of the graph.  

The function `energy_edges_fragment` implements the edge removal process, and takes as argument the graph to remove, and a Dict which `fragmentation` entry represents the amount of energy to remove from the graph.

The function `get_fragmentation_probability_energy_edge_removal` computes the probability of fragmentation given an amount of energy to remove by using a Monte Carlo method: the previously described procedure is repeated to approximate the fragmentation probability.
The number of iterations and the edge removal probability are given as parameters.

The function `get_fragmentation_energy_threshold_edge` approximates the energy percolation threshold (i.e., the “energy” to remove in order to obtain a probability of fragmentation of 0.5) by using the same bisection method as `get_fragmentation_probability_threshold_node`, using `get_fragmentation_probability_energy_edge_removal` for each step.

### Removing nodes
The “energy of a node” is defined as the “energy” needed to remove it, i.e as the sum of the energies of all edges adjacent to it. The method for removing nodes is similar to the one for removing edges: the probability weight of each node is proportional to the inverse of its energy. Here the energy of each node is stored as a networkx attribute. The main difference being that as a node gets removed, the energies of the neighbouring nodes have to be decreased by the energy of the edges that were linked to them via the removed node. This process can also leave isolated nodes that have an energy of zero, with undefined probability weights. This situation is avoided by removing isolated neighbours as well as the chosen node.  
Note that the `get_fragmentation_energy_threshold_node` and `get_fragmentation_probability_energy_node_removal` functions are, respectively, similar to `get_fragmentation_energy_threshold_edge` and `get_fragmentation_probability_energy_edge_removal` in terms of parameters.

The function `energy_nodes_fragment` implements the node removal process and has the same argments as the `energy_edges_fragment` function.

Note that for these functions to work properly, the `energy` attribute of the nodes needs to be initialized before passing the graph to the functions. The `init_nodes_energy` function sets the nodes attribute as previously described, given a graph where the edge energy attributes are already defined. See the examples for more details.

## Hole size detection
The `capsidgraph.generator` modules provides methods to compute the statistic destribution of "hole sizes" in graph under fragmentation.
### Definition of hole size
Let $G = (V,E)$ be the graph we want to fragment, let $G' = (V',E')$ the graph where nodes or edges have been removed ($G \neq G'$ and $V \neq \emptyset$), we want to define the size of the largest hole in $G'$.

Consider the set of connected components of maximal size $\Set{C_0=(V_0, E_0),...,C_p=(V_p,E_p)}$

 let $i \in \Set{0,...,p}$ and $\bar{C_i}=(\bar{V_i}, \bar{E_i})$ where

 $\bar{V_i} = V \setminus V_i$
 
 $\bar{E_i} = \Set{\Set{u,v} :\Set{u,v}\in E, u \in \bar{V_i}, v \in \bar{V_i}}$

let $H_i$ be the size of> the largest connected component of $C_i$

Then the largest hole size of $G'$, $H(G') = \min_{0 \leq j \leq p }{H_j}$

We say that a graph is _hole-fragmented_ if $H(G') \geq \lfloor \frac{|V|}{2} \rfloor$

If $G=G'$, by convention we define $H(G')=0$
Similarly, if $V'=\emptyset$ we define $H(G')=|V|$

### Implementations
The `get_hole_size` function implements the above-mentioned algorithm to determine the size of the largest hole of a graph.

The `get_hole_size_distribution` function estimate the probability distribution of the "Hole size" random variable. The fragment function used to remove edges or nodes from the graph is passed as a parameter.



# Graph generation
The `capsidgraph.generator` module provides tools for generating the graphs representing the protein subunit interaction networks of viral capsids. This module can generate graphs representing the interaction network of a viral capsid.

## Icosahedral capsid generation
### Generation of a triangular face
A face in the capsid graph is generated by starting with a list of edges forming a tile, as well as the two translation vectors $\vec{T_x}$ and $\vec{T_y}$ that translate this tile such that the translated copies tile the plane with regularly spaced 6-fold symmetry axis ((0,0) being one of them). A rescaled equilateral triangular face is then "cut out" from this tiling as follows: The three vertices of the triangle must be chosen to coincide with 6-fold symmetry axes. The position of this triangle is defined by two integers $(h,k)$ such that $(h \vec{T_x}, k \vec{T_y})$ is the vector from (0,0) to another vertex of the triangle. The function `create_icosahedral_face_edges` automates this algorithm. The three parameters required are the edges defining the tile, the translation vectors $\vec{T_x}$ and $\vec{T_y}$, and the values for $h$ and $k$.  
For instance, the following code generated a face tiled with triangles.  
```python
    from capsidgraph.generator import create_icosahedral_face_edges
    tile = [
        ((0,0),(0,1)),
        ((0,0),(-1,1)),
        ((0,0),(-1,0)),
        ((0,0),(0,-1)),
        ((0,0),(1,-1)),
        ((0,0),(1,0)),
    ]
    Tx = (1, 0)
    Ty = (0,1)
    face_edges, axis = create_icosahedral_face_edges(tile,Tx,Ty,1,2)
```
The tile drawn in Cartesian coordinates can be seen below  
<img src="img/graph1.png" height=300 alt="The edges of a single tile">  
Which gives the following tiling  
<img src="img/graph2.png" height=300 alt="The corresponding tiling">  
The function returns only the edges inside the defined triangle  
<img src="img/graph3.png" height=300 alt="The edges of a triangular face">  
which corresponds to the following edges and vertices delimiting the face
```python
faceEdges = [((1, 0), (1, 1)), ((1, 1), (1, 2)), ((1, -1), (1, 0)), ((2, 0), (2, 1)), ((2, -1), (2, 0)), ((1, 0), (0, 1)), ((1, 1), (0, 2)), ((2, 0), (1, 1)), ((2, -1), (1, 0)), ((3, -1), (2, 0)), ((1, 0), (0, 0)), ((1, 1), (0, 1)), ((2, 0), (1, 0)), ((2, 1), (1, 1)), ((3, 0), (2, 0)), ((1, 0), (1, -1)), ((1, 1), (1, 0)), ((1, 2), (1, 1)), ((2, 0), (2, -1)), ((2, 1), (2, 0)), ((0, 1), (1, 0)), ((0, 2), (1, 1)), ((1, 0), (2, -1)), ((1, 1), (2, 0)), ((2, 0), (3, -1)), ((0, 0), (1, 0)), ((0, 1), (1, 1)), ((1, 0), (2, 0)), ((1, 1), (2, 1)), ((2, 0), (3, 0))]
faceVertices = ((0, 0), (1, 2), (3, -1))
```

Some patterns are predefined in the `generator.icosahedral_patterns` module.

### Graph generation
Given one triangular face of an icosahedron as described above, the graph of an icosahedral surface lattice is obtained as follows. By copying this face 20 times and merging nodes such that those 20 triangles match at the icosahedral edges, we generate the graph of the corresponding capsid. This construction is done with the help of a dictionary that associates to each node the face(s) associated with it, along with its coordinates inside each face. Initially every node only belongs to one triangular face, and this information is stored in terms of its coordinates inside that face. To "glue" two triangular faces $F_0$ and $F_1$ together, we first rotate all the points of $F_0$ by 60 degrees along one of the triangular vertices. Then for each point in $F_0$, we search for a point in $F_1$ with the same coordinates. If such a point is found, those points are merged. For this, they are retained as a single node in the graph, and their coordinate lists in the dictionary are concatenated.  

<img src="img/faceMerge.png" height=300 alt="Two faces being fused together.">

Once all 20 triangular faces have been correctly assembled into an icosahedral surface lattice, the associated graph is fully assembled. In our example, the graph looks like this (graph visualised with the `draw` function of networkx)  

<img src="img/graph4.png" height=300 alt="The final graph">  
This construction is implemented in the function `create_icosahedral_capsid_graph` which takes the following parameters as input:
*  `face_edges` : a list of edges corresponding to the edges of a triangular face;
* `triangle_vertices` : the three vertices delimiting the faces (tuple of 3 points);
* `bond_strength` : optional, this list needs to have the same length as `face_edges`. This parameter defines the bond energy of every vertex of the capsid graph, where  `bond_strength[i]` corresponds to the energy of the bond `face_edges[i]`.

## Cubic graph generation
Some nanoparticle interaction networks can be constructed through a process similar to the icosahedral generation algorithm, where the icosahedral pattern is replaced by a cubic one. The faces are now squares which vertices sit in 4-fold symetry axis of the lattice.
The `create_cubic_face_edges` function is similar to `create_icosahedral_face_edges`, with the exeption of the h and k parameters that have been replaced with the coordinates of one vertex of the square face (given that one of them is sitting at the origin this defined the square).

Some patterns are predefined in the `generator.cubic_patterns` module.

### Example
The following code
```python
from capsidgraph.generator import (
    cubic_patterns,
    create_cubic_face_edges,
    create_cubic_capsid_graph,
)
P = cubic_patterns.AALS_48_PATTERN
[edges, Tx, Ty, face_side_edge] = P
face_edges, face_square_vertices = create_cubic_face_edges(edges, Tx, Ty, face_side_edge)
G = create_cubic_capsid_graph(face_edges, face_square_vertices)
```
will create a graph G, isomorphic to the following graph :

<img src="img/aalas48.png" height="400" style="margin:auto">
