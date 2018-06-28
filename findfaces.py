#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 15:22:51 2018

@author: eion + michelle
"""

import networkx as nx
from itertools import chain
import numpy as np
import math
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

def find_faces(G, nbr_lists):
    edges = ((u,v) for u,v,d in G.edges(data=True))
    verse = ((v,u) for u,v,d in G.edges(data=True))
    darts = list(chain(edges, verse))
#    print(len(darts))
    faces = set()
    used = set()
    for start in darts:
        if start not in used:
            face = []
            hunt = True
            dart = start
            while hunt:
                used.add(dart)
                face.append(dart[0])
                if dart == start and face != [dart[0]]:
                    hunt = False
                else:
                    pivot = dart[1]
                    print(pivot)
                    nbr_list = nbr_lists[pivot]
                    index = nbr_list.index(dart[0])
                    new_index = (index+1)%len(nbr_list)
                    dart = (pivot, nbr_list[new_index])
            faces.add(tuple(set(face)))
#    print(len(faces))
    return faces



G=nx.read_gexf('northcarolina.gexf')

nbr_lists=genclockwiseneighbors(G)
w = csv.writer(open("nc_clockwisenbrs.csv","w"))
for key, val in nbr_lists.items():
    w.writerow([key, val])

faces = find_faces(G, nbr_lists)