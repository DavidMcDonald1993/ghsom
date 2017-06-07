
# coding: utf-8

# In[3]:

import numpy as np
import networkx as nx

from Queue import Queue
from threading import Thread
from threading import current_thread

from sklearn.manifold import MDS
from sklearn.metrics.pairwise import euclidean_distances

def distances(X, d):
    return euclidean_distances(X[:,:d])

def stress(distances, Ds, tri):
    return np.sqrt(np.sum((distances[tri] - np.sqrt(Ds[tri])) ** 2) / np.sum(Ds[tri]))
    
def classic_mds(Ds):
        
    #centering matrix
    n = len(Ds)
    C = np.identity(n) - np.ones((n, n)) / n
    
    #double centered similarity matrix
    K = - 0.5 * np.dot(np.dot(C, Ds), C)
    
    #eigen decompose K
    l, U = np.linalg.eigh(K)
    
    mask = l > 0
    l = l[mask]
    U = U[:, mask]
    
    ##sort eigenvalues (and eigenvectors) into ascending order
    idx = l.argsort()[::-1]
    l = l[idx]
    U = U[:,idx]
    
#     if k < 1:
        
#         print "attempting to preserve 95% of variance"
        
#         total = sum(l)
#         k = 0
#         var = sum(l[:k]) / total
#         while var < variance_preserved:
#             k += 1
#             var = sum(l[:k]) / total

#     print "embedding to {} dimensions".format(k)
    
    #position matrix
    X = np.dot(U, np.diag(l ** 0.5))
    
    #upper triangle
    tri = np.triu_indices(n=len(Ds), k = 1, m=len(Ds))

    #determine the k that minimised stress
    stresses = np.array([stress(distances(X, d), Ds, tri) for d in range(1, min(50, U.shape[1]))])
    k = stresses.argmin() + 1
    
    print "embedding to {} dimensions with stress {}".format(k, stresses[k-1])
    
    return X[:,:k], stresses

##function to generate benchmark graph
def benchmark_hierarchical_graph(edge_path, c1_path, c2_path):

    #construct graph from edge list
    G = nx.read_edgelist(edge_path)

    #create dictionarys of attributes
    c1 = read_attributes(c1_path)
    c2 = read_attributes(c2_path)

    #set attributes of G
    nx.set_node_attributes(G, 'firstlevelcommunity', c1)
    nx.set_node_attributes(G, 'secondlevelcommunity', c2)
    
    #return graph
    return G

##function to generate benchmark graph
def benchmark_graph(edge_path, c_path):
    
    #construct graph from edge list
    G = nx.read_edgelist(edge_path)

    #create dictionarys of attributes
    c = read_attributes(c_path)

    #set attributes of G
    nx.set_node_attributes(G, 'firstlevelcommunity', c)
    
    #return graph
    return G

##function to read in attributes from file and return a dictionary
def read_attributes(filepath):
    
    #initialise dictionary
    d = {}
    
    #open file
    with open(filepath) as f:
        
        #iterate over lines in f
        for l in f:
            
            #separate into key and value
            k, v = l.split()
            
            #add to dictionary
            d[k] = v
    
    #return
    return d

def filter_nodes_with_no_embedding(G, D):
    
    for n in G.nodes():
    
        if n not in D.keys():
            print "{} not in D, removing it from the network".format(n)
            G.remove_node(n)

##save embedding to graph
def set_embedding(G, X):
    
    embedding = {k: v for k, v in zip(G.nodes(), X)}
    
    nx.set_node_attributes(G, "embedding", embedding)
    
    
## get embedding
def get_embedding(G):

    return np.array([v for k, v in nx.get_node_attributes(G, "embedding").items()])
    
def main_hierarchical(network, first_level, second_level, filename):
    
    #import graph from file
    G = benchmark_hierarchical_graph(network, first_level, second_level)
    
    #only embed largest subgraph
    H = max(nx.connected_component_subgraphs(G), key=len)
    
    #embed into X
    D = nx.floyd_warshall(H)
    
    #remove nodes from H with no embedding
    filter_nodes_with_no_embedding(H, D)    

    ##to array
    Ds = np.array([[D[i][j] ** 2 for j in D.keys()] for i in D.keys()])

    #mds
    X, stresses = classic_mds(Ds)
    
    #save embedding to nodes of G
    set_embedding(H, X)
    
    #write gml file
    nx.write_gpickle(H, filename)

    return X, Ds, stresses
    
def main_binary(network, first_level, filename):
    
    #import graph from file
    G = benchmark_graph(network, first_level)
    
    #only embed largest subgraph
    H = max(nx.connected_component_subgraphs(G), key=len)
    
    #embed into X
    D = nx.floyd_warshall(H)
    
    #remove nodes from H with no embedding
    filter_nodes_with_no_embedding(H, D)   
    
    ##to array
    Ds = np.array([[D[i][j] ** 2 for j in D.keys()] for i in D.keys()])
    
    #mds
    X, stresses = custom_mds(Ds) 

    
#     print "embedding to {} dimensions".format(k)

#     mds = MDS(n_components = k, dissimilarity = "precomputed", n_jobs=-1, max_iter=10000)
#     X = mds.fit_transform(D)
    
    #save embedding to nodes of G
    set_embedding(H, X)
    
    #write gml file
    nx.write_gpickle(H, filename)
    
    return X, Ds, stresses

def main(txt, filename, D=None, delimiter="\t"):
    
    #import graph from file
    G = nx.read_edgelist(txt, delimiter=delimiter)
    
    #only embed largest subgraph
    H = max(nx.connected_component_subgraphs(G), key=len)
    
    #embed into X
    if D == None:
        #if no precomuped distance matrix then use floyd
        D = nx.floyd_warshall(H)
    
    #remove nodes from H with no embedding
    filter_nodes_with_no_embedding(H, D) 
    
    ##array of squared shortest path distances
    Ds = np.array([[D[i][j] ** 2 for j in D.keys()] for i in D.keys()])
    
    print "computed distance matrix"
    
#     #validation testing for best dimension
#     stresses = np.array([])
    
#     for k in range(min_k, max_k + 1):
        
#         mds = MDS(n_components = k, dissimilarity = "precomputed", n_jobs=-1, max_iter=10000)
#         mds.fit(D)

#         stress = mds.stress_
#         X = mds.embedding_
        
#         print stress

#         stresses = np.append(stresses, stress)
        
#         if len(stresses) > 1:
#             gradient = np.gradient(stresses, 5)[-1]

#             print "k={}, stress={}, gradient={}".format(k, stress, gradient)

#             if np.abs(gradient) < eps:
#                 print "BREAK"
#                 break
    
#     mds = MDS(n_components = k, dissimilarity = "precomputed", n_jobs=-1, max_iter=10000)
#     mds.fit(D)

#     stress = mds.stress_
#     X = mds.embedding_
    X, stresses = classic_mds(Ds)
    
    #save embedding to nodes of G
    set_embedding(H, X)

    #write gml file
    nx.write_gpickle(H, filename)

    return X, Ds, stresses


# In[3]:

X, D, stresses = main("Uetz_screen.txt", "embedded_yeast_uetz.gpickle", delimiter="\t")


# In[4]:

X, D, stresses = main_hierarchical("benchmarks/network.dat", "benchmarks/community_first_level.dat", 
                  "benchmarks/community_second_level.dat", "benchmarks/hierarchical_benchmark.gpickle")


# In[ ]:

stresses


# In[34]:

G = nx.read_gpickle("benchmarks/hierarchical_benchmark.gpickle")


# In[35]:

X = get_embedding(G)


# In[36]:

X.shape


# In[37]:

G.nodes(data=True)


# In[ ]:



