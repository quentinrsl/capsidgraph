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

#Affiche l'image d'une face icosahédrique de la capsid
def createFaceImage(pattern:Tuple[List[Edge],Point,Point,int],h:int,k:int,lineWidth:float=10,lineColor:Tuple[int,int,int]=(0,0,0),scale:float=500):

    #On récupère un pattern avec les vecteurs de translations pour le pavage
    [edges, Tx, Ty] = pattern[:3]
    
    #ON récupère les coordonées du traingle dans lequel on va découper l'image 
    A,B,C = createTriangle(h,k,Tx,Ty)
    theta= atan2(carth(B)[0],carth(B)[1])   #Angle duquel on va faire tourner toute l'image pour l'avoir avec la bonne orientation
    #Largeur et hauteur de l'image finale
    width = int(rotate(carth(C),theta)[0]) * scale
    height = int(rotate(carth(B),theta)[1]) * scale
    #On duplique le pattern autant que nécessaire
    for i in range(h+k):
        edges = extend(edges,Tx)
    for i in range(h+k):
        edges = extend(edges,Ty)
    #Liste des segments en coordonnées hexagonales
    edges = extractTriangle(edges,h,k,Tx,Ty)
    
    #On les converti en coordonnées carthésiennes, et on les tourne de theta
    carthEdges = []
    for P1,P2 in edges:
        carthEdges.append((rotate(carth(P1),theta), rotate(carth(P2),theta)))
    
    #On les dessine sur un canvas
    im = Image.new(mode="RGBA", size=(width,height), color=(255,255,255))
    draw = ImageDraw.Draw(im)
    for (x1,y1),(x2,y2) in carthEdges:
        draw.line([x1 * scale,y1 * scale,x2 * scale,y2 * scale], fill=lineColor,width = lineWidth)
    
    draw.polygon([(0,0),(width, height/2),(width,0)],(255,255,255))
    draw.polygon([(0,height),(width, height/2),(width,height)],(255,255,255))
    im.show()
