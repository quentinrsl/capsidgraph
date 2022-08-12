from graph_generator import *
from math import sqrt,cos,sin,atan2
from PIL import Image, ImageDraw

Point = Tuple[float,float]
Edge = Tuple[Point,Point]

#Convert hexagonal coordiantes into cathesian ones
def carth(P:Point)->Point:
    (x,y) = P
    return (
        (x+y/2), 
        (y*sqrt(3)/2)
    )

#Rreturn point rotated of theta around the origin
def rotate(P:Point,theta:float)->Point:
    (x,y) = P
    return (
        x*cos(theta) - y*sin(theta), 
        x*sin(theta) + y*cos(theta)
    )

#return the image of an icosahedral face of the capsid
def getFaceImage(pattern:Tuple[List[Edge],Point,Point,int],h:int,k:int,lineWidth:float=10,lineColor:Tuple[int,int,int]=(0,0,0),scale:float=500)->Image:

    #Get the required data to generate the tiling
    [edges, Tx, Ty] = pattern[:3]
    
    #Coordinates of the triangle to cut into the tiling
    A,B,C = createTriangle(h,k,Tx,Ty)
    theta= atan2(carth(B)[0],carth(B)[1])   #Angle by which we need to rotate the image to get the right orientation
    #width and height of the final image
    width = int(rotate(carth(C),theta)[0]) * scale
    height = int(rotate(carth(B),theta)[1]) * scale
    #Duplicate the tile as much as necessary
    for i in range(h+k):
        edges = extend(edges,Tx)
    for i in range(h+k):
        edges = extend(edges,Ty)
    #Get the list of extyracted lines
    edges = extractTriangle(edges,h,k,Tx,Ty)
    
    #Convert them back into carthesian coordinates and rotate them
    carthEdges = []
    for P1,P2 in edges:
        carthEdges.append((rotate(carth(P1),theta), rotate(carth(P2),theta)))
    
    #Draw the lines and show the image
    im = Image.new(mode="RGBA", size=(width,height), color=(255,255,255))
    draw = ImageDraw.Draw(im)
    for (x1,y1),(x2,y2) in carthEdges:
        draw.line([x1 * scale,y1 * scale,x2 * scale,y2 * scale], fill=lineColor,width = lineWidth)
    
    draw.polygon([(0,0),(width, height/2),(width,0)],(255,255,255))
    draw.polygon([(0,height),(width, height/2),(width,height)],(255,255,255))
    return im
