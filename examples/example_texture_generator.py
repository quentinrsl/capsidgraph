from capsidgraph.generator import icosahedral_patterns, create_icosahedral_texture


from PIL import Image
P = icosahedral_patterns.PATTERN_6434
img = create_icosahedral_texture(P, 2, 1)
img2 = Image.open("tests/test_texture.png")