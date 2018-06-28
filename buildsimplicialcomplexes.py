"""
Created on Thu Jun 28 15:22:51 2018

@author: eion + michelle + austin + corey
"""

import networkx as nx
from itertools import chain
import numpy as np
import math
import csv


def find_angle(y,x):
    """Find angle of a vector with respect to x axis"""
    theta = math.atan(y/x)
    if x<0:
        theta = theta + math.pi
    return theta


def gen_cclockwise_neighbors(G,xstring,ystring):
    """Puts neighbors of every vertex in clockwise order

    :param G: networkx graph with nodes indexed by GeoID, as well as coordinates for each node
    :param xstring: string key for x coordinate in G
    :param ystring: string key for y coordinate in G
    :return: dictionary where keys are GeoIDs for each node, and values are a list of neighbors in
        counterclockwise order
    """
    cclockwisenbrlist={}
    for node in G.nodes():
        coord = np.array([float(G.nodes[node][xstring]),float(G.nodes[node][ystring])])
        neighbors = G.neighbors(node)
        angles = []
        nbrs = []
        for n in neighbors:
            nbrcoord = np.array([float(G.nodes[n][xstring]), float(G.nodes[n][ystring])])
            u = nbrcoord-coord
            theta = find_angle(u[1],u[0])
            angles.append(theta)
            nbrs.append(str(n))
        cclockwisenbrlist[node]=[n for _,n in sorted(zip(angles,nbrs))]
    return cclockwisenbrlist


def find_faces(G, nbr_lists):
    """Finds faces of a mostly planar graph G

    :param G: networkx graph
    :param nbr_lists: dictionary mapping GeoIDs to lists of neighbors in counter-clockwise order
    :return: dictionary mapping indexed faces to a list of GeoIDs for vertices of that face in counter-clockwise order
    """
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
    """For a face with more than 3 vertices, breaks face into triangles

    :param face: string ID for the face
    :param verts: list of GeoIDs of vertices in the face in counterclockwise order
    :return: trianglelist, a dictionary mapping triangle IDs to a list of vertices in each triangle, where triangles
        roughly triangulate the face
    """
    trianglelist={}
    numverts = len(verts)
    if numverts > 3:
        for subidx in range(numverts-2):
            trianglelist[face+str(subidx)]=[verts[0],verts[subidx+1],verts[subidx+2]]
    else:
        trianglelist[face]=verts
    return trianglelist

def get_edge_dict(G):
    """assign each edge in a graph a string index and get the GeoIDs of endpoints

    :param G: networkx graph indexed by GeoID
    :return: dictionary mapping unique edge indexes to the GeoIDs of their endpoints
    """
    edge_dict={}
    idx = 0
    for e in list(G.edges()):
        verts = list(e)
        edge_dict["E"+str(idx)]=verts
        idx+=1
    return edge_dict

def get_dims_and_idx(G,edgedict,faces):
    """Build dictionaries that give an order for the entry of every simplex into the filtered simplicial complex

    :param G: networkx graph
    :param edgedict: dictionary mapping unique edge indexes to the GeoIDs of their endpoints
    :param faces: dictionary mapping indexed faces to a list of GeoIDs for vertices of that face in
        counter-clockwise order
    :return:
        dims: dictionary mapping the ID of each simplex to its dimension
        indexes: dictionary mapping the ID of each simplex to the order in which it enters the simplex
        tridict: dictionary mapping triangles (2-simplexes) to the list of GeoIDs of their vertices
    """
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
    """from two vertices, find the edge ID of the edge between them

    :param verts: list of two GeoIDs
    :param edgedict: dictionary mapping edge IDs to GeoIDs of their endpoints
    :return: edge ID of the edge between the vertices given by verts
    """
    for e, everts in edgedict.items():
        if set(everts) == set(verts):
            return e
    print('Error: No Match Found.')
    return 0


def write_simplicial_complex(filename,dims,indexes, edgedict, tridict):
    """write filtered simplicial complex to .dat file

    :param filename: filename of the output file
    :param dims: dictionary mapping simplex IDs to their dimensions
    :param indexes: dictionary giving order in which simplex IDs enter the filtered simplicial complex
    :param edgedict: dictionary mapping edge IDs to GeoIDs of their endpoints
    :param tridict: dictionary mapping triangles (2-simplexes) to the list of GeoIDs of their vertices
    :return: 
    """
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
nbr_list = gen_cclockwise_neighbors(G,"INTPTLAT","INTPTLON")
edict = get_edge_dict(G)
faces = find_faces(G,nbr_list)
dims,idxs,tridict = get_dims_and_idx(G,edict,faces)
write_simplicial_complex("wytestsc.dat",dims,idxs,edict,tridict)
