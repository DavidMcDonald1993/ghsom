
# coding: utf-8

# In[1]:

import numpy as np
import networkx as nx
import som_functions as som
import math
import sklearn.metrics as met
import sklearn.manifold as man
from time import time

##floyds embedding
def floyd_embedding(G):
    
    n = len(G)
    
    fl = nx.floyd_warshall(G)
    
    #intitialise distance matrix
    D = np.zeros((n, n))
    
    #find closest k neighbours
    for i in range(n):
        n1 = G.nodes()[i]
        for j in range(n):
            n2 = G.nodes()[j]
            D[i,j] = fl[n1][n2]
    
    C = np.identity(n) - np.ones((n, n)) / n
    
    #similarity matrix
    K = - 1/2 * np.dot(np.dot(C, D ** 2), C)
    
    #eigen decompose K
    l, U = np.linalg.eigh(K)
    
    ##sort eigenvalues (and eigenvectors) into ascending order
    idx = l.argsort()[::-1]
    l = l[idx]
    U = U[:,idx]
    
    s = sum(l)
    
    k = len(l)
    var = 1
    
    while var > 0.95:
        k -= 1
        var = sum(l[:k]) / s
    
    k += 1
    
    #position matrix
    X = np.dot(U[:,:k], np.diag(l[:k] ** 0.5))
    
    return X

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

##save embedding to graph
def set_embedding(G, X):
    
    #get number of niodes in the graph
    num_nodes = nx.number_of_nodes(G)
    
    #dimension of embedding
    d = len(X[0])
    
    #iterate over a;; the nodes and save their embedding
    for i in range(num_nodes):
        for j in range(d):
            G.node[G.nodes()[i]]['embedding'+str(j)] = X[i,j]    

def main_hierarchical(network, first_level, second_level):
    
    #import graph from file
    G = benchmark_hierarchical_graph(network, first_level, second_level)
    
    #only embed largest subgraph
    H = max(nx.connected_component_subgraphs(G), key=len)
    
    #embed into X
    X = floyd_embedding(H)
    
    #save embedding to nodes of G
    set_embedding(H, X)
    
    #write gml file
    nx.write_gml(H, 'embedded_network_{}.gml'.format(network.split('_')[0]))

def main_binary(network, first_level, gml_filename):
    
    #import graph from file
    G = benchmark_graph(network, first_level)
    
    #only embed largest subgraph
    H = max(nx.connected_component_subgraphs(G), key=len)
    
    #embed into X
    X = floyd_embedding(H)
    
    #save embedding to nodes of G
    set_embedding(H, X)
    
    #write gml file
    nx.write_gml(H, gml_filename)

def main(txt, gml_filename):
    
    #import graph from file
    G = nx.read_edgelist(txt)
    
    #only embed largest subgraph
    H = max(nx.connected_component_subgraphs(G), key=len)
    
    #embed into X
    X = floyd_embedding(H)
    
    #save embedding to nodes of G
    set_embedding(H, X)
    
    #write gml file
    nx.write_gml(H, gml_filename)


# In[4]:

G = nx.read_gml("polbooks.gml")

#only embed largest subgraph
H = max(nx.connected_component_subgraphs(G), key=len)

#embed into X
X = floyd_embedding(H)

#save embedding to nodes of G
set_embedding(H, X)

#write gml file
nx.write_gml(H, "embedded_pollbooks.gml")


# In[ ]:



