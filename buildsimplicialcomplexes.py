"""
Created on Thu Jun 28 15:22:51 2018

@author: eion + michelle + austin + cory
"""

import networkx as nx
from itertools import chain
import numpy as np
import math
import csv
import dionysus as d
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.lines as mlines
import matplotlib.cm
import numpy as np
from scipy.spatial.distance import pdist


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

def gettriangles(face,verts,edgedict):
    """For a face with more than 3 vertices, breaks face into triangles

    :param face: string ID for the face
    :param verts: list of GeoIDs of vertices in the face in counterclockwise order
    :param edgedict: dictionary mapping unique edge indexes to the GeoIDs of their endpoints
    :return: trianglelist, a dictionary mapping triangle IDs to a list of vertices in each triangle, where triangles
        roughly triangulate the face
    """
    trianglelist={}
    numverts = len(verts)
    newedges = []
    if numverts > 3 and numverts < 40:
        for subidx in range(numverts-2):
            # print(verts[0],verts[subidx+1],verts[subidx+2])
            if(len(set([verts[0],verts[subidx+1],verts[subidx+2]])))==3:
            # if True:
                trianglelist[face+"T"+str(subidx)]=[verts[0],verts[subidx+1],verts[subidx+2]]
                if lookup_edge_key([verts[0],verts[subidx+2]],edgedict)==0:
                    edgedict["E"+face+"T"+str(subidx)]=[verts[0],verts[subidx+2]]
                    # print([verts[0],verts[subidx+2]])
                    newedges.append("E"+face+"T"+str(subidx))
    else:
        if len(set(verts))==3:
        # if True:
            trianglelist[face]=verts
    return trianglelist,newedges

def triangulate_faces(facedict):
    triangledict = {}
    idx = 0
    for f in facedict:
        if len(set(facedict[f])) == 3:
            triangledict["T" + str(idx)] = facedict[f]
            idx += 1
        elif len(set(facedict[f])) > 3:
            verts = facedict[f]
            numverts = len(verts)
            for subidx in range(numverts-2):
                testtri = [verts[0], verts[subidx+1], verts[subidx+2]]
                if len(set(testtri)) == 3:
                    triangledict["T" + str(idx)] = facedict[f]
                    idx += 1
    return triangledict

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
    done_nodes = []
    tridict={}
    time={}
    for currdist in range(maxdist+1):
        nodes = [key for key,value in disttobd.items() if value==currdist]
        for node in nodes:
            dims[node] = 0
            time[node] = currdist
            indexes[node] = idx
            done_nodes.append(node)
            idx+=1
            for e in list(rem_edges):
                everts = rem_edges[e]
                if len(set(everts) & set(done_nodes))==2:
                    dims[e]=1
                    indexes[e] = idx
                    time[e] = currdist
                    rem_edges.pop(e)
                    idx+=1
            for f in list(rem_faces):
                fverts = rem_faces[f]
                if len(set(fverts) & set(done_nodes)) == len(set(fverts)):
                    trilist, newedges = gettriangles(f,fverts,edgedict)
                    for e in newedges:
                        dims[e]=1
                        indexes[e] = idx
                        time[e] = currdist
                        idx+=1
                    for t, tverts in trilist.items():
                        dims[t] = 2
                        indexes[t] = idx
                        time[t] = currdist
                        idx+=1
                        tridict[t] = tverts
                    rem_faces.pop(f)
    return dims, indexes, tridict, time

def get_dims_and_idx_from_multiple(Gedge, faces, maxfiltration, numfiltrations):
    filtrationwidth = maxfiltration/numfiltrations
    nodes = Gedge.nodes()
    rem_edges = nx.get_edge_attributes(Gedge,'weight')
    rem_faces = faces.copy()
    done_edges = {}
    idx = 0
    dims = {}
    idxs = {}
    times = {}
    for n in nodes:
        dims[n] = 0
        idxs[n] = idx
        times[n] = 0
        idx += 1
    for currdist in range(numfiltrations):
        currfilt = (currdist+1)*filtrationwidth
        for e in list(rem_edges.keys()):
            if rem_edges[e] < currfilt:
                dims[e] = 1
                idxs[e] = idx
                times[e] = currdist
                idx += 1
                done_edges[e]=rem_edges[e]
                rem_edges.pop(e)
        for f in list(rem_faces.keys()):
            if check_face(rem_faces[f], done_edges) == True:
                # if e == ('11001009601','11001009602'):
                #     print("Here")
                fidx = tuple(rem_faces[f])
                dims[fidx] = 2
                idxs[fidx] = idx
                times[fidx] = currdist
                idx += 1
                rem_faces.pop(f)
    return dims, idxs, times

def lookup_edge_key_2(e, edgedict):
    for edge in edgedict.keys():
        if set(e) == set(edge):
            return edge
    return 0

def lookup_edge_key(verts,edgedict):
    """from two vertices, find the edge ID of the edge between them

    :param verts: list of two GeoIDs
    :param edgedict: dictionary mapping edge IDs to GeoIDs of their endpoints
    :return: edge ID of the edge between the vertices given by verts
    """
    for e, everts in edgedict.items():
        if set(everts) == set(verts):
            return e
    return 0

def check_face(faceverts,edgedict):
    if lookup_edge_key_2([faceverts[0], faceverts[1]], edgedict) != 0 and lookup_edge_key_2([faceverts[1], faceverts[2]], edgedict) != 0 and lookup_edge_key_2([faceverts[0], faceverts[2]], edgedict) != 0:
        return True
    return False

def build_filtered_complex(dims, indexes, edgedict, tridict, time):
    sortedkeys = sorted(indexes, key=indexes.get)
    f = d.Filtration()
    for key in sortedkeys:
        if dims[key]==0:
            f.append(d.Simplex([indexes[key]],time[key]))
        elif dims[key]==1:
            verts = edgedict[key]
            vertices = [indexes[v] for v in verts]
            f.append(d.Simplex(vertices, time[key]))
        elif dims[key]==2:
            verts = tridict[key]
            vertices = [indexes[v] for v in verts]
            f.append(d.Simplex(vertices, time[key]))
    return f

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

def plot_simplexes(G,xstring,ystring,edgedict,triangledict):
    """generate plot of colorcoded simplexes by distance from boundary

    :param G: network which must contain distance to boundary
    :param xstring:
    :param ystring:
    :param edgedict: dictionary of all edges and contained vertices
    :param triangledict: dictionary of all triangles and contained vertices
    :return:
    """
    colorlist = []
    for n in G.nodes():
        G.node[n]['pos'] = (float(G.node[n][xstring]), float(G.node[n][ystring]))
        colorlist.append(float(G.node[n]['dist_to_boundary']))

    pos = nx.get_node_attributes(G, 'pos')
    nppos = np.array([val for k, val in pos.items()])
    xmin, ymin = np.min(nppos, 0)
    xmax, ymax = np.max(nppos, 0)
    disttobd = nx.get_node_attributes(G, 'dist_to_boundary')
    maxdist = max(list(disttobd.values()))

    clist = ['k'] * (maxdist + 1)
    clist[0] = 'r'
    clist[1] = 'b'
    clist[2] = 'g'
    clist[3] = 'c'
    clist[4] = 'm'
    clist[5] = 'y'
    fig, ax = plt.subplots()
    for currdist in range(maxdist + 1):
        boundary = [k for k, val in time.items() if val == currdist]
        boundarytri = [k for k in boundary if dims[k] == 2]
        boundaryedges = [k for k in boundary if dims[k] == 1]
        boundaryverts = [k for k in boundary if dims[k] == 0]
        patches = []
        points = np.array([np.array(pos[k]) for k in boundaryverts])
        plt.scatter(points[:, 0], points[:, 1], s=5, c=clist[currdist], alpha=0.4)
        for e in boundaryedges:
            verts = edgedict[e]
            beg = np.array(pos[verts[0]])
            end = np.array(pos[verts[1]])
            line = mlines.Line2D(np.array([beg[0], end[0]]), np.array([beg[1], end[1]]), lw=1., alpha=0.4,
                                 c=clist[currdist])
            ax.add_line(line)
        for t in boundarytri:
            poly = np.zeros((3, 2))
            verts = triangledict[t]
            for i in np.arange(3):
                poly[i, :] = np.array(pos[verts[i]])
            triangle = Polygon(poly, True)
            patches.append(triangle)
        p = PatchCollection(patches, alpha=0.4)
        p.set_color(clist[currdist])
        ax.add_collection(p)
    plt.axis([xmin, xmax, ymin, ymax])

    plt.show(block=True)

def plot_simplexes_from_multiple(G_edge,dims,times):
    """generate plot of colorcoded simplexes by distance from boundary

    :param G: network which must contain distance to boundary
    :param xstring:
    :param ystring:
    :param edgedict: dictionary of all edges and contained vertices
    :param triangledict: dictionary of all triangles and contained vertices
    :return:
    """

    pos = nx.get_node_attributes(G_edge, 'pos')
    nppos = np.array([val for k, val in pos.items()])
    xmin, ymin = np.min(nppos, 0)
    xmax, ymax = np.max(nppos, 0)
    maxdist = max(list(times.values()))

    cmap = matplotlib.cm.get_cmap('viridis')
    clist = [cmap(val) for val in np.linspace(0,1,maxdist+1)]
    fig, ax = plt.subplots()
    for currtime in range(maxdist):
        boundary = [k for k, val in times.items() if val == currtime]
        boundarytri = [k for k in boundary if dims[k] == 2]
        boundaryedges = [k for k in boundary if dims[k] == 1]
        boundaryverts = [k for k in boundary if dims[k] == 0]
        patches = []
        if boundaryverts:
            points = np.array([np.array(pos[k]) for k in boundaryverts])
            plt.scatter(points[:, 0], points[:, 1], s=5, alpha=0.4, c=clist[currtime])
        for e in boundaryedges:
            verts = list(e)
            beg = np.array(pos[verts[0]])
            end = np.array(pos[verts[1]])
            line = mlines.Line2D(np.array([beg[0], end[0]]), np.array([beg[1], end[1]]), lw=1., alpha=0.4,
                                 c=clist[currtime])
            ax.add_line(line)
        for t in boundarytri:
            poly = np.zeros((3, 2))
            verts = list(t)
            for i in np.arange(3):
                poly[i, :] = np.array(pos[verts[i]])
            triangle = Polygon(poly, True)
            patches.append(triangle)
        p = PatchCollection(patches, alpha=0.4)
        p.set_color(clist[currtime])
        ax.add_collection(p)
    plt.axis([xmin, xmax, ymin, ymax])

    plt.show(block=True)

def compute_ph_boundary_using_graph_distance(gexffile,districtid,latstring,lonstring):
    """ Computes the persistent homology of an boundary outward march with graph distance

    :param gexffile: location of input gexf
    :param districtid: identifier for the district whose boundary we're examining
    :param latstring: string ID for latitude field of centroid
    :param lonstring: string ID for longitude field of centroid
    :return:
    """
    G = nx.read_gexf(gexfile)

    get_dist_to_boundary(G,districtid)

    nbr_list = gen_cclockwise_neighbors(G,lonstring,latstring)
    edict = get_edge_dict(G)
    faces = find_faces(G,nbr_list)
    dims,idxs,tridict, time = get_dims_and_idx(G,edict,faces)


    filtration = build_filtered_complex(dims, idxs, edict, tridict, time)

    plot_simplexes(G,lonstring,latstring,edict,tridict)

    m = d.homology_persistence(filtration)

    dgms = d.init_diagrams(m,filtration)
    d.plot.plot_bars(dgms[1], show=True)

def build_adj_matrix_percent(datafile, keystring, xstring, ystring, attr_list):
    data = np.zeros((0,6))
    G = nx.Graph()
    idxs = {}
    idx = 0
    with open(datafile) as geoid:
        reader = csv.reader(geoid,delimiter=',')
        headers = next(reader)
        key_idx = headers.index(keystring)
        x_idx = headers.index(xstring)
        y_idx = headers.index(ystring)
        dat_idx = [headers.index(dstring) for dstring in attr_list]
        for row in reader:
            race = np.array([[float(row[i]) for i in dat_idx]])
            race = race/np.sum(race)
            data = np.concatenate((data,race))
            G.add_node(row[key_idx], pos=(float(row[x_idx]),float(row[y_idx])))
            idxs[idx] = row[key_idx]
            idx += 1
    condenseddist = pdist(data)
    num_nodes = len(G.nodes())
    for i in range(num_nodes-1):
        for j in range(num_nodes-i-1):
            condensedidx = int(num_nodes*j - j*(j+1)/2 + i - 1 - j)
            G.add_edge(idxs[i], idxs[i+j+1], weight = condenseddist[condensedidx])
    return G



Gspatial = nx.read_gexf('gexf/DC_tracts.gexf')
Gdem = build_adj_matrix_percent('rawdata/DC_2010Census_Race.csv','GEOID','xcoord','ycoord',
                                ['Tract_2010Census_DP1_DP0090001',
                                 'Tract_2010Census_DP1_DP0090002',
                                 'Tract_2010Census_DP1_DP0090003',
                                 'Tract_2010Census_DP1_DP0090004',
                                 'Tract_2010Census_DP1_DP0090005',
                                 'Tract_2010Census_DP1_DP0090006'])
pos = nx.get_node_attributes(Gdem, 'pos')
xcoord={k:u for k,(u,v) in pos.items()}
ycoord={k:v for k,(u,v) in pos.items()}

nx.set_node_attributes(Gspatial, pos, 'pos')
nx.set_node_attributes(Gspatial, xcoord, 'xcoord')
nx.set_node_attributes(Gspatial, ycoord, 'ycoord')
nbr_list = gen_cclockwise_neighbors(Gspatial,'xcoord', 'ycoord')
faces_spatial = find_faces(Gspatial,nbr_list)
tris_spatial = triangulate_faces(faces_spatial)
edict_dem = get_edge_dict(Gdem)

dims, idxs, times = get_dims_and_idx_from_multiple(Gdem, tris_spatial, .5, 20)
plot_simplexes_from_multiple(Gdem, dims, times)
# nx.draw(Gspatial, pos=nx.get_node_attributes(Gspatial,'pos'),node_size=10)
# nx.draw(Gdem, pos=nx.get_node_attributes(Gdem,'pos'),node_size=10)
# plt.show(block=True)
# get_dims_and_idx_from_multiple(G, edict, eweight, faces, 1., 1.)
