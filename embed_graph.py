
# coding: utf-8

# In[18]:

import numpy as np
import networkx as nx

from Queue import Queue
from threading import Thread
from threading import current_thread

from sklearn.manifold import MDS
    
# def custom_mds(D, k=5, variance_preserved=0.95):
        
#     #centering matrix
#     n = len(D)
#     C = np.identity(n) - np.ones((n, n)) / n
    
#     #similarity matrix
#     K = - 0.5 * np.dot(np.dot(C, D ** 2), C)
    
#     #eigen decompose K
#     l, U = np.linalg.eigh(K)
    
#     ##sort eigenvalues (and eigenvectors) into ascending order
#     idx = l.argsort()[::-1]
#     l = l[idx]
#     U = U[:,idx]
    
#     if k < 1:
    
#         s = (l - l[-1]) / (l[0] - l[-1])
#         print l
#         print s

#         total = sum(s)
#         k = 0
#         var = sum(s[:k]) / total
#         while var < variance_preserved:
#             k += 1
#             var = sum(s[:k]) / total
        
#         print "attempting to preserve 95% of variance".format(k)
    
#     print "embedding to {} dimensions".format(k)
    
    
#     #position matrix
#     X = np.dot(U[:,:k], np.diag(l[:k] ** 0.5))
    
#     return X

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
    D = np.array([[D[i][j] for j in D.keys()] for i in D.keys()])

    k = 10
    
    print "embedding to {} dimensions".format(k)

    mds = MDS(n_components = k, dissimilarity = "precomputed", n_jobs=-1, max_iter=1000)
    X = mds.fit_transform(D)
    
    #save embedding to nodes of G
    set_embedding(H, X)
    
    #write gml file
    nx.write_gpickle(H, filename)

def main_binary(network, first_level, filename, k = 5):
    
    #import graph from file
    G = benchmark_graph(network, first_level)
    
    #only embed largest subgraph
    H = max(nx.connected_component_subgraphs(G), key=len)
    
    #embed into X
    D = nx.floyd_warshall(H)
    
    #remove nodes from H with no embedding
    filter_nodes_with_no_embedding(H, D)   
    
    ##to array
    D = np.array([[D[i][j] for j in D.keys()] for i in D.keys()])
    
    #mds
    X = custom_mds(D, k = k) 

    
#     print "embedding to {} dimensions".format(k)

#     mds = MDS(n_components = k, dissimilarity = "precomputed", n_jobs=-1, max_iter=10000)
#     X = mds.fit_transform(D)
    
    #save embedding to nodes of G
    set_embedding(H, X)
    
    #write gml file
    nx.write_gpickle(H, filename)

def main(txt, filename, D=None, delimiter="\t", min_k = 1, max_k = 50, eps = 10, k=10):
    
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
    
    ##to array
    D = np.array([[D[i][j] for j in D.keys()] for i in D.keys()])
    
    print "computed distance matrix"
    
#     #validation testing for best dimension
#     stresses = np.array([])
    
#     for k in [10, 20, 100, 250, 500, 1000]:
        
#         mds = MDS(n_components = k, dissimilarity = "precomputed", n_jobs=-1, max_iter=10000)
#         mds.fit(D)

#         stress = mds.stress_
#         X = mds.embedding_

#         stresses = np.append(stresses, stress)
        
#         print "k={}, stress={}".format(k, stress)
        
#         if len(stresses) > 1:
#             gradient = np.gradient(stresses, 5)[-1]

#             print "k={}, stress={}, gradient={}".format(k, stress, gradient)

#             if np.abs(gradient) < eps:
#                 print "BREAK"
#                 break
    
    mds = MDS(n_components = k, dissimilarity = "precomputed", n_jobs=-1, max_iter=10000)
    mds.fit(D)

    stress = mds.stress_
    X = mds.embedding_
    
    
    #save embedding to nodes of G
    set_embedding(H, X)

    #write gml file
    nx.write_gpickle(H, filename)
    
    return k, stress, X


# In[7]:

G = nx.read_edgelist("reactome_edgelist.txt")
H = max(nx.connected_component_subgraphs(G), key=len)
D = nx.floyd_warshall(H)


# In[19]:

stresses = main("reactome_edgelist.txt", "embedded_yeast_reactome.gpickle", D = D, k = 20)


# In[ ]:

stresses


# In[22]:

import matplotlib.pyplot as plt


# In[24]:

plt.plot(stresses)
plt.show()


# In[4]:

main_binary("benchmarks/network.dat", "benchmarks/community.dat", "embedded_benchmark.gpickle", k = 10)


# In[9]:

G = nx.read_gpickle("embedded_yeast_reactome.gpickle")


# In[10]:

G.nodes(data=True)


# In[ ]:



