import networkx as nx
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import csv

def findangle(y,x):
    theta = math.atan(y/x)
    if x<0:
        theta = theta + math.pi
    return theta

def genclockwiseneighbors(G):
    clockwisenbrlist={}
    for node in G.nodes():
        coord = np.array([float(G.nodes[node]['INTPTLAT10']),float(G.nodes[node]['INTPTLON10'])])
        neighbors = G.neighbors(node)
        angles = []
        nbrs = []
        for n in neighbors:
            nbrcoord = np.array([float(G.nodes[n]['INTPTLAT10']), float(G.nodes[n]['INTPTLON10'])])
            u = nbrcoord-coord
            theta = findangle(u[1],u[0])
            angles.append(theta)
            nbrs.append(str(n))
        clockwisenbrlist[node]=[n for _,n in sorted(zip(angles,nbrs))]
    return clockwisenbrlist

G=nx.read_gexf('northcarolina.gexf')

clockwisenbrlist=genclockwiseneighbors(G)
w = csv.writer(open("nc_clockwisenbrs.csv","w"))
for key, val in clockwisenbrlist.items():
    w.writerow([key, val])