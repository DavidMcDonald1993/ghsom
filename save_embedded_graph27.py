
# coding: utf-8

# In[2]:

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


# In[8]:

# #import graph from file
# # G = benchmark_hierarchical_graph('1_network.dat', '1_community_first_level.dat', '1_community_second_level.dat')
# G = nx.karate_club_graph()

# #only embed largest subgraph
# H = max(nx.connected_component_subgraphs(G), key=len)

# #embed into X
# X = floyd_embedding(H)


# In[9]:

# %matplotlib notebook

# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(X[:,0], X[:,1], X[:,2])
# plt.show()


# In[3]:

# # generate graph
# # # G = nx.karate_club_graph()
# G = nx.read_gml('football.gml')
# # # G = nx.read_gml('adjnoun.gml')
# G = nx.read_gml('dolphins_labelled.gml')
# # # G = nx.read_gml('polbooks.gml')
# # # G = benchmark_graph('bin_network.dat', 'community.dat')
# # # G = benchmark_hierarchical_graph('network.dat', 
# # #                                  'community_first_level.dat', 
# # #                                  'community_second_level.dat')
# # # G = benchmark_hierarchical_graph('literature_network.dat', 
# # #                                  'literature_community_first_level.dat', 
# # #                                  'literature_community_second_level.dat')
# # # G = benchmark_hierarchical_graph('literature_network_32.dat', 
# # #                                  'literature_community_first_level_32.dat', 
# # #                                  'literature_community_second_level_32.dat')
# # # G = benchmark_hierarchical_graph('literature_network_double.dat', 
# # #                                  'literature_community_first_level_double.dat', 
# # #                                  'literature_community_second_level_double.dat')
# # # G = benchmark_hierarchical_graph('rand_network.dat', 
# # #                                  'rand_community_first_level.dat', 
# # #                                  'rand_community_second_level.dat')
# # # G = benchmark_hierarchical_graph('1000network.dat', 
# # #                                  '1000community_first_level.dat', 
# # #                                  '1000community_second_level.dat')
# # G = benchmark_hierarchical_graph('yang_network.dat', 
# #                                  'yang_community_first_level.dat', 
# #                                  'yang_community_second_level.dat')

# # print('loaded G',len(G),'nodes')

# labels = ['club']
# # # labels = ['value']
# labels = ['group']
# # # labels = ['firstlevelcommunity']
# # labels = ['firstlevelcommunity','secondlevelcommunity']

# # #divide graph into connected components
# connected_graphs = connected_components(G)
# print('calculated number of connected components')

# # #embedding
# # X = dsd_embedding(H)
# X = floyd_embedding(G)
    
# # #save embedding to graph
# set_embedding(G, X)
# print('embedded graph')

# # #save graph
# nx.write_gml(G, "embedded_dolphin.gml")
    
# print('saved graph')


# In[24]:

# import networkx as nx
# import numpy as np
# import os

# def graph_measures(G):
    
#     return nx.degree_assortativity_coefficient(G), nx.density(G), nx.node_connectivity(G)

# os.chdir('C:\Miniconda3\Jupyter\GHSOM_simplex_dsd')

# num_nodes = 64

# G = nx.read_gml('embedded_karate.gml')
# print graph_measures(G)
# G = nx.read_gml('embedded_dolphin.gml')
# print graph_measures(G)
# G = nx.read_gml('embedded_polbooks.gml')
# print graph_measures(G)
# G = nx.read_gml('embedded_football.gml')
# print graph_measures(G)
# print

# for i in range(10):
#     G = nx.barabasi_albert_graph(num_nodes, np.random.randint(1, 16))

#     print graph_measures(G)
# #     print nx.degree_pearson_correlation_coefficient(G)


# In[25]:

# os.chdir('C:\Miniconda3\Jupyter\GHSOM_simplex_dsd\parameter_tests_small_density_simulated_annealing')

# for i in range(200):
    
#     G = nx.read_gml('embedded_network_{}.gml'.format(i))
    
#     print graph_measures(G)


# In[ ]:



