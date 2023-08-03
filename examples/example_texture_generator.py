from capsidgraph.generator import icosahedral_patterns, create_icosahedral_texture
from PIL import Image

"""
This example shows how to use capsidgraph.generator to create an icosahedral texture that would fit on a face of an icosahedron from a lattice parttern.
"""

P = icosahedral_patterns.PATTERN_6434
img = create_icosahedral_texture(P, 2, 1)
img.show()