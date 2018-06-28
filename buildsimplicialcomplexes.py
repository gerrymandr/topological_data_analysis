"""
Created on Thu Jun 28 15:22:51 2018

@author: eion + michelle + austin + corey
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


def genclockwiseneighbors(G,latstring,lonstring):
    clockwisenbrlist={}
    for node in G.nodes():
        coord = np.array([float(G.nodes[node][latstring]),float(G.nodes[node][lonstring])])
        neighbors = G.neighbors(node)
        angles = []
        nbrs = []
        for n in neighbors:
            nbrcoord = np.array([float(G.nodes[n][latstring]), float(G.nodes[n][lonstring])])
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
#     faces = set()
    faces = {}
    faceidx = 0
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
                    nbr_list = nbr_lists[pivot]
                    index = nbr_list.index(dart[0])
                    new_index = (index+1)%len(nbr_list)
                    dart = (pivot, nbr_list[new_index])
            # faces.add(tuple(face[:len(face)-1]))
            faces["F"+str(faceidx)]=face[:len(face)-1]
            faceidx+=1
#    print(len(faces))
    return faces


def get_boundary_nodes(G, district):
    #takes in VTD adjacency graph G and district identifier (string),
    #returns list of boundary nodes in that district
    complement_nodes = []
    for node in G.nodes():
        if G.node[node]['DISTRICT'] != district:
            complement_nodes.append(node)
    boundary_nodes_of_district = nx.node_boundary(G, complement_nodes)
    return boundary_nodes_of_district

def get_dist_to_boundary(G, district):
    #computes dist to boundary of vertex_set (list of vertices) for each node in G,
    #adds this dist as a new attribute (dist_to_boundary) to each node in G
    nx.set_node_attributes(G,len(G.edges())+1,'dist_to_boundary')
    boundary_nodes_of_district = get_boundary_nodes(G, district)
    for source_node in boundary_nodes_of_district:
        for target_node in G.nodes():
            current_dist = G.node[target_node]['dist_to_boundary']
            updated_dist = nx.shortest_path_length(G,source_node,target_node)
            G.node[target_node]['dist_to_boundary'] = updated_dist if updated_dist < current_dist else current_dist
    return G

def gettriangles(face,verts):
    trianglelist={}
    numverts = len(verts)
    if numverts > 3:
        for subidx in range(numverts-2):
            trianglelist[face+str(subidx)]=[verts[0],verts[subidx+1],verts[subidx+2]]
    else:
        trianglelist[face]=verts
    return trianglelist

def get_edge_dict(G):
    edge_dict={}
    idx = 0
    for e in list(G.edges()):
        verts = list(e)
        edge_dict["E"+str(idx)]=verts
        idx+=1
    return edge_dict

def get_dims_and_idx(G,edgedict,faces):
    disttobd = nx.get_node_attributes(G,'dist_to_boundary')
    maxdist = max(list(disttobd.values()))
    idx = 0
    dims = {}
    indexes = {}
    rem_edges=edgedict.copy()
    rem_faces=faces.copy()
    nodes = list(G)
    done_nodes = []
    tridict={}
    for currdist in range(maxdist):
        nodes = [key for key,value in disttobd.items() if value==currdist]
        for node in nodes:
            dims[node] = 0
            indexes[node] = idx
            nodes.remove(node)
            done_nodes.append(node)
            idx+=1
            for e in list(rem_edges):
                everts = rem_edges[e]
                if len(set(everts) & set(done_nodes))==2:
                    dims[e]=1
                    indexes[e] = idx
                    rem_edges.pop(e)
                    idx+=1
            for f in list(rem_faces):
                fverts = rem_faces[f]
                if len(set(fverts) & set(done_nodes)) == len(fverts):
                    trilist = gettriangles(f,fverts)
                    for t, tverts in trilist.items():
                        dims[t] = 2
                        indexes[t] = idx
                        idx+=1
                        tridict[t] = tverts
                    rem_faces.pop(f)
    return dims, indexes, tridict

def lookup_edge_key(verts,edgedict):
    for e, everts in edgedict.items():
        if set(everts) == set(verts):
            return e
    print('Error: No Match Found.')
    return 0


def write_simplicial_complex(filename,dims,indexes, edgedict, tridict):
    sortedkeys = sorted(indexes, key=indexes.get)
    file = open(filename,"w")
    for key in sortedkeys:
        if dims[key]==0:
            file.write("0\n")
        elif dims[key]==1:
            verts = edgedict[key]
            file.write("1 "+str(indexes[verts[0]])+" "+str(indexes[verts[1]])+"\n")
        elif dims[key]==2:
            verts = tridict[key]
            edges = []
            edges.append(lookup_edge_key(verts[0:2],edgedict))
            edges.append(lookup_edge_key(verts[1:3],edgedict))
            edges.append(lookup_edge_key([verts[0],verts[2]],edgedict))
            file.write("2 "+str(indexes[edges[0]])+" "+str(indexes[edges[1]])+" "+str(indexes[edges[2]])+"\n")

    file.close()

G = nx.read_gexf('wytestgraph2.gexf')
get_dist_to_boundary(G,'00')
nbr_list = genclockwiseneighbors(G,"INTPTLAT","INTPTLON")
edict = get_edge_dict(G)
faces = find_faces(G,nbr_list)
dims,idxs,tridict = get_dims_and_idx(G,edict,faces)
write_simplicial_complex("testsc.dat",dims,idxs,edict,tridict)
