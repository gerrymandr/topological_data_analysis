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
    if x!=0:
        theta = math.atan(y/x)
    elif y>0:
        theta = math.pi/2
    else:
        theta = math.pi*3/2
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

def gettriangles(face,verts,G):
    """For a face with more than 3 vertices, breaks face into triangles

    :param face: string ID for the face
    :param verts: list of GeoIDs of vertices in the face in counterclockwise order
    :param G: networkx graph
    :return:
        trianglelist, a dictionary mapping triangle IDs to a list of vertices in each triangle, where triangles
        roughly triangulate the face
        newedges, list of edge IDs that have been added to the graph for triangulation reasons
    """
    trianglelist={}
    numverts = len(verts)
    newedges = []
    if numverts > 3 and numverts < 40:
        for subidx in range(numverts-2):
            if(len(set([verts[0],verts[subidx+1],verts[subidx+2]])))==3:
                trianglelist[face+"T"+str(subidx)]=[verts[0],verts[subidx+1],verts[subidx+2]]
                if lookup_edge_key([verts[0],verts[subidx+2]],G)==0:
                    G.add_edge(verts[0],verts[subidx+2])
                    e = lookup_edge_key([verts[0],verts[subidx+2]],G)
                    newedges.append(e)
    else:
        if len(set(verts))==3:
            trianglelist[face]=verts
    return trianglelist,newedges

def triangulate_faces(facedict):
    '''triangulate all faces which are not triangular

    :param facedict: dictionary listing unique face IDs with vertices involved in them
    :return: dictionary of unique triangle IDs with their vertices
    '''
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
                    triangledict["T" + str(idx)] = testtri
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

def get_dims_and_idx(G,edgelist,faces):
    """Build dictionaries that give an order for the entry of every simplex into the filtered simplicial complex

    :param G: networkx graph
    :param edgedict: dictionary mapping unique edge indexes to the GeoIDs of their endpoints
    :param faces: dictionary mapping indexed faces to a list of GeoIDs for vertices of that face in
        counter-clockwise order
    :return:
        dims: dictionary mapping the ID of each simplex to its dimension
        indexes: dictionary mapping the ID of each simplex to the order in which it enters the simplex
        time: dictionary mapping ID of each simplex to time it enters the simplex
    """
    disttobd = nx.get_node_attributes(G,'dist_to_boundary')
    maxdist = max(list(disttobd.values()))
    idx = 0
    dims = {}
    indexes = {}
    rem_edges=list(edgelist)
    rem_faces=faces.copy()
    done_nodes = []
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
                everts = list(e)
                if len(set(everts) & set(done_nodes))==2:
                    dims[e]=1
                    indexes[e] = idx
                    time[e] = currdist
                    rem_edges.remove(e)
                    idx+=1
            for f in list(rem_faces):
                fverts = rem_faces[f]
                if len(set(fverts) & set(done_nodes)) == len(set(fverts)):
                    trilist, newedges = gettriangles(f,fverts,G)
                    for e in newedges:
                        dims[e]=1
                        indexes[e] = idx
                        time[e] = currdist
                        idx+=1
                    for t, tverts in trilist.items():
                        tkey=tuple(tverts)
                        dims[tkey] = 2
                        indexes[tkey] = idx
                        time[tkey] = currdist
                        idx+=1
                    rem_faces.pop(f)
    return dims, indexes, time

def get_dims_and_idx_from_multiple(Gedge, faces, maxfiltration, numfiltrations):
    '''Build dictionaries for dimensions, entry time with edges from one graph and faces from another

    :param Gedge: graph with shared vertices and the edges
    :param faces: list of faces from another graph with the vertices of each face
    :param maxfiltration: maximum epsilon for building epsilon balls in the simplicial complex
    :param numfiltrations: number of filtration steps
    :return:
        dims: dimensions of each simplex
        idx: order that each simplex enters the simplicial complex
        times: timestep at which each simplex enters the simplicial complex
    '''
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
                fidx = tuple(rem_faces[f])
                dims[fidx] = 2
                idxs[fidx] = idx
                times[fidx] = currdist
                idx += 1
                rem_faces.pop(f)
    return dims, idxs, times

def lookup_edge_key_dict(e, edgedict):
    """lookup whether an edge exists in a dictionary of edges with keys given by tuple of vertices

    :param e: list of vertices in the potential edge
    :param edgedict: dictionary of edges (can have any value)
    :return:
        edge: edge id of the found edge
    """
    for edge in edgedict.keys():
        if set(e) == set(edge):
            return edge
    return 0

def lookup_edge_key(verts,G):
    """from two vertices, find the edge ID of the edge between them

    :param verts: list of two GeoIDs
    :param G: networkx object with all edges
    :return: edge ID of the edge between the vertices given by verts
    """
    for edge in G.edges():
        if set(edge) == set(verts):
            return edge
    return 0

def check_face(faceverts,edgedict):
    """check whether all the edges of a face are an a given list of edges

    :param faceverts: vertices of the face in order
    :param edgedict: dictionary containing all permissible edges
    :return:
        True if all edges are contained, False otherwise
    """
    if lookup_edge_key_dict([faceverts[0], faceverts[1]], edgedict) != 0 \
            and lookup_edge_key_dict([faceverts[1], faceverts[2]], edgedict) != 0 \
            and lookup_edge_key_dict([faceverts[0], faceverts[2]], edgedict) != 0:
        return True
    return False

def build_filtered_complex(dims, indexes, time):
    """Build a filtered complex from dictionaries giving the dimensions, subsimplices, and times of entry

    :param dims: dictionary with keys giving simplex, values giving dimension
    :param indexes: Dictionary with keys giving the simplex, value giving the order in which simplex enters
    :param time: Dictionary with keys giving simplex, value giving time of entry
    :return:
    """
    sortedkeys = sorted(indexes, key=indexes.get)
    f = d.Filtration()
    for key in sortedkeys:
        if dims[key]==0:
            f.append(d.Simplex([indexes[key]],time[key]))
        elif dims[key]==1:
            verts = list(key)
            vertices = [indexes[v] for v in verts]
            f.append(d.Simplex(vertices, time[key]))
        elif dims[key]==2:
            verts = list(key)
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

def plot_simplexes_from_boundary(G, xstring, ystring, time, dims):
    """generate plot of colorcoded simplexes by distance from boundary

    :param G: network which must contain distance to boundary
    :param xstring:
    :param ystring:
    :param dims: dictionary of simplex -> dimension
    :param times: dictionary of simplex -> entry time
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

    cmap = matplotlib.cm.get_cmap('viridis')
    clist = [cmap(val) for val in np.linspace(0, 1, maxdist + 1)]
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
            verts = list(e)
            beg = np.array(pos[verts[0]])
            end = np.array(pos[verts[1]])
            line = mlines.Line2D(np.array([beg[0], end[0]]), np.array([beg[1], end[1]]), lw=1., alpha=0.4,
                                 c=clist[currdist])
            ax.add_line(line)
        for t in boundarytri:
            poly = np.zeros((3, 2))
            verts = list(t)
            for i in np.arange(3):
                poly[i, :] = np.array(pos[verts[i]])
            triangle = Polygon(poly, True)
            patches.append(triangle)
        p = PatchCollection(patches, alpha=0.4)
        p.set_color(clist[currdist])
        ax.add_collection(p)
    plt.axis([xmin, xmax, ymin, ymax])
    plt.interactive(False)
    plt.show()

def plot_simplexes_from_multiple(G_edge,dims,times):
    """generate plot of colorcoded simplexes by distance from boundary

    :param G: network which must contain distance to boundary
    :param xstring:
    :param ystring:
    :param dims: dictionary of simplex -> dimension
    :param times: dictionary of simplex -> entry time
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
    for currtime in range(maxdist+1):
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
    plt.interactive(False)
    plt.show()

def compute_ph_boundary_using_graph_distance(gexfile,districtid,latstring,lonstring):
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
    faces = find_faces(G,nbr_list)
    dims, idxs, time = get_dims_and_idx(G,G.edges(),faces)

    filtration = build_filtered_complex(dims, idxs, time)

    plot_simplexes_from_boundary(G, lonstring, latstring, time, dims)

    m = d.homology_persistence(filtration)

    dgms = d.init_diagrams(m,filtration)
    d.plot.plot_bars(dgms[1], show=True)

def build_adj_matrix_percent(datafile, keystring, xstring, ystring, attr_list):
    """Read in a csv with demographic percentile data

    :param datafile: location of csv
    :param keystring: GEOID
    :param xstring:
    :param ystring:
    :param attr_list: list of the field IDs for the attributes being examined
    :return: graph with demographic data attached and edge weights from pairwise distance
    """
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

def compute_ph_dem_spatial_from_percent(spatial_gexf, dem_csv, idstring, xstring, ystring, attrstrings, numfils=20):
    """runs persistent homology on the simplicial complex with edges from one graph and faces from another. Requires
    demographic data to be percentile

    :param spatial_gexf: path to gexf of rook adjacency
    :param dem_csv: path to csv with demographic data
    :param idstring:
    :param xstring:
    :param ystring:
    :param attrstrings: list of desired demographic attributes
    :param numfils:
    :return:
    """
    Gspatial = nx.read_gexf(spatial_gexf)
    Gdem = build_adj_matrix_percent(dem_csv, idstring, xstring, ystring, attrstrings)
    pos = nx.get_node_attributes(Gdem, 'pos')
    xcoord={k:u for k,(u,v) in pos.items()}
    ycoord={k:v for k,(u,v) in pos.items()}

    nx.set_node_attributes(Gspatial, pos, 'pos')
    nx.set_node_attributes(Gspatial, xcoord, 'xcoord')
    nx.set_node_attributes(Gspatial, ycoord, 'ycoord')
    nbr_list = gen_cclockwise_neighbors(Gspatial,'xcoord', 'ycoord')
    faces_spatial = find_faces(Gspatial,nbr_list)
    tris_spatial = triangulate_faces(faces_spatial)
    weights = nx.get_edge_attributes(Gdem,'weight')
    maxfil = np.median([w for w in weights.values()])/10.

    dims, idxs, times = get_dims_and_idx_from_multiple(Gdem, tris_spatial, maxfil, numfils)
    plot_simplexes_from_multiple(Gdem, dims, times)

    filtration = build_filtered_complex(dims, idxs, times)

    m = d.homology_persistence(filtration)

    dgms = d.init_diagrams(m,filtration)
    d.plot.plot_bars(dgms[1], show=True)

def run_simple_example():
    """runs a simple illustrative example"""
    Gspatial = nx.Graph()
    Gspatial.add_nodes_from(range(9))
    Gspatial.add_edges_from(
        [(0, 1), (0, 2), (1, 3), (2, 3), (3, 4), (3, 6), (3, 8), (4, 5), (5, 6), (5, 7), (6, 7), (7, 8)])
    xcoord = {0: 0., 1: 1., 2: 0., 3: 1., 4: 4., 5: 3., 6: 2., 7: 2., 8: 1., }
    ycoord = {0: 0., 1: 0., 2: -1., 3: -1., 4: -1., 5: -2., 6: -2., 7: -3., 8: -4.}
    pos = {k: np.array([xcoord[k], ycoord[k]]) for k in xcoord.keys()}
    nx.set_node_attributes(Gspatial, xcoord, 'xcoord')
    nx.set_node_attributes(Gspatial, ycoord, 'ycoord')
    nx.set_node_attributes(Gspatial, pos, 'pos')

    Gdem = nx.Graph()
    Gdem.add_nodes_from(range(9))
    Gdem.add_weighted_edges_from(
        [(0, 1, .1), (0, 2, .1), (0, 3, .2), (1, 3, .1), (2, 3, 0.1), (3, 4, .1), (4, 8, 0.1), (3, 8, .1), (5, 6, .2),
         (5, 7, .1), (6, 7, .1)])
    nx.set_node_attributes(Gdem, xcoord, 'xcoord')
    nx.set_node_attributes(Gdem, ycoord, 'ycoord')
    nx.set_node_attributes(Gdem, pos, 'pos')
    nbr_list = gen_cclockwise_neighbors(Gspatial, 'xcoord', 'ycoord')
    nbr_list = {k: [int(i) for i in v] for k, v in nbr_list.items()}
    faces_spatial = find_faces(Gspatial, nbr_list)
    tris_spatial = triangulate_faces(faces_spatial)
    weights = nx.get_edge_attributes(Gdem, 'weight')

    dims, idxs, times = get_dims_and_idx_from_multiple(Gdem, tris_spatial, 1., 20)
    plot_simplexes_from_multiple(Gdem, dims, times)

    filtration = build_filtered_complex(dims, idxs, times)

    m = d.homology_persistence(filtration)

    dgms = d.init_diagrams(m,filtration)
    d.plot.plot_bars(dgms[1], show=True)

### EXAMPLE FOR COMPUTING PERSISTENT HOMOLOGY USING A DISTANCE FROM BOUNDARY FILTRATION
# compute_ph_boundary_using_graph_distance('gexf/wytestgraph2.gexf','00','INTPTLAT','INTPTLON')

### EXAMPLE FOR COMPUTING THE PERSISTENT HOMOLOGY OF A SIMPLICIAL COMPLEX WITH EDGES GIVEN BY DEMOGRAPHIC DISTANCE AND
### SIMPLICES GIVEN BY SPATIAL DISTANCE
# compute_ph_dem_spatial_from_percent('gexf/DC_tracts.gexf',
#                                     'rawdata/DC_2010Census_Race.csv',
#                                     'GEOID',
#                                     'xcoord',
#                                     'ycoord',
#                                     ['Tract_2010Census_DP1_DP0090001',
#                                      'Tract_2010Census_DP1_DP0090002',
#                                      'Tract_2010Census_DP1_DP0090003',
#                                      'Tract_2010Census_DP1_DP0090004',
#                                      'Tract_2010Census_DP1_DP0090005',
#                                      'Tract_2010Census_DP1_DP0090006'])

### EXTREMELY SIMPLE EXAMPLE
# run_simple_example()