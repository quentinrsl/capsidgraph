from capsidlib.graph_generator import createCapsidGraph


#Create the capsid of the polyomavirus with given energy on a, b and c edges
def createPolyomavirusCapsid(Ea,Eb,Ec):
    faceEdges = [((1,0),(1,1)), #c
        ((1,0),(1,-1)),#b
        ((1,0),(0,1)),#b
        ((1,0),(2,-1)),#c
        ((1,0),(0,0)),#a
        ((1,0),(2,0)),#c
        ((1,1),(2,0)),#b
        ((1,-1),(2,-1)),#c
        ((2,0),(2,1)),#a
        ((2,0),(2,-1)),#c
        ((2,0),(3,-1)),#c
        ((2,0),(3,0)),#b
        ((2,-1),(2,-2)),#b
        ((2,-1),(3,-2)),#a
        ((2,-1),(3,-1))#b
    ]
    axis = ((0, 0), (2, 1), (3, -2))

    energy =   [Ec,Eb,Eb,Ec,Ea,Ec,Eb,Ec,Ea,Ec,Ec,Eb,Eb,Ea,Eb]
    return createCapsidGraph(faceEdges,axis,bondStrength=energy)
