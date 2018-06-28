import networkx as nx

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
    for currdist in range(maxdist):
        nodes = [key for key,value in disttobd.items() if value==currdist]
        for node in nodes:
            dims[node] = 0
            indexes[node] = idx
            nodes.remove(node)
            done_nodes.append(node)
            idx+=1
            for e, everts in rem_edges.items():
                if len(set(everts) & set(done_nodes))==2:
                    dims[e]=1
                    indexes[e] = idx
                    rem_edges.pop(e)
                    idx+=1
            for f, fverts in rem_faces.items():
                if len(set(fverts) & set(done_nodes)) == len(fverts):
                    trilist = gettriangles(f,fverts)
                    for t, tverts in trilist:
                        dims[t] = 2
                        indexes[t] = idx
                        idx+=1
                    rem_faces.pop(f)
    return dims, indexes

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
            file.write("0")
        else if dims[key]==1:
            verts = edgedict[key]
            file.write("1 "+str(indexes[verts[0]])+" "+str(index[verts[1]]))
        else if dims[key]==2:
            verts = tridict[key]
            edges = []
            edges.append(lookup_edge_key(verts[0:2],edgedict))
            edges.append(lookup_edge_key(verts[1:3],edgedict))
            edges.append(lookup_edge_key([verts[0],verts[2]],edgedict))
            file.write("2 "+str(indexes[edges[0]])+" "+str(indexes[edges[1]])+" "+str(indexes[edges[2]]))

    file.close()

