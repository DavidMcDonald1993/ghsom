
# coding: utf-8

# In[1]:

import numpy as np
import networkx as nx
import som_functions as som
import math
import sklearn.metrics as met
import sklearn.manifold as man
from time import time

# ##floyds embedding
# def floyd_embedding(G):
    
#     n = len(G)
    
#     fl = nx.floyd_warshall(G)
    
#     #intitialise distance matrix
#     D = np.zeros((n, n))
    
#     #find closest k neighbours
#     for i in range(n):
#         n1 = G.nodes()[i]
#         for j in range(n):
#             n2 = G.nodes()[j]
#             D[i,j] = fl[n1][n2]
            
#     return D
    
def custom_mds(distance_dict):
    
    #construct distance matrix D from dictionary
    D = np.array([[distance_dict[i][j] for j in distance_dict[i]] for i in distance_dict])
        
    #centering matrix
    n = len(distance_dict)
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
    
    #link nodes to embedding
    X = {k: v for k, v in zip(distance_dict, X)}
    
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

def filter_nodes_with_no_embedding(G, D):
    
    for n in G.nodes():
    
        if n not in D:
            print "{} not in D, removing it from the network".format(n)
            G.remove_node(n)

##save embedding to graph
def set_embedding(G, X):
    
    #get number of niodes in the graph
#     num_nodes = nx.number_of_nodes(G)
    
    #dimension of embedding
    d = len(X[G.nodes()[0]])
    
    #iterate over a;; the nodes and save their embedding
    for n in G.nodes():
        for j in range(d):
            G.node[n]['embedding{}'.format(j)] = X[n][j]    

def main_hierarchical(network, first_level, second_level):
    
    #import graph from file
    G = benchmark_hierarchical_graph(network, first_level, second_level)
    
    #only embed largest subgraph
    H = max(nx.connected_component_subgraphs(G), key=len)
    
    #embed into X
    D = nx.floyd_warshall(H)
    
    #remove nodes from H with no embedding
    filter_nodes_with_no_embedding(H, D)    
    
    #mds
    X = custom_mds(D) 
    
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
    D = nx.floyd_warshall(H)
    
    #remove nodes from H with no embedding
    filter_nodes_with_no_embedding(H, D)    
    
    #mds
    X = custom_mds(D) 
    
    #save embedding to nodes of G
    set_embedding(H, X)
    
    #write gml file
    nx.write_gml(H, gml_filename)

def main(txt, gml_filename, D=None):
    
    #import graph from file
    G = nx.read_edgelist(txt)
    
    #only embed largest subgraph
    H = max(nx.connected_component_subgraphs(G), key=len)
    
    #embed into X
    if D == None:
        #if no precomuped distance matrix then use floyd
        D = nx.floyd_warshall(H)
    
    #remove nodes from H with no embedding
    filter_nodes_with_no_embedding(H, D)    
    
    #mds
    X = custom_mds(D)   
    
    #save embedding to nodes of G
    set_embedding(H, X)
    
    #write gml file
    nx.write_gml(H, gml_filename)


# In[72]:

main_binary("benchmarks/network.dat", "benchmarks/community.dat", "benchmarks/embedded_benchmark_3.gml")


# In[3]:

graph_file = "Y2H_union.txt"
gml_filename = "embedded_yeast_union_rel.gml"
distance_file = "yeast_union_rel_similarity_GOSim.csv"

labels = np.genfromtxt(distance_file, delimiter=',', usecols=0, dtype=str)
raw_data = np.genfromtxt(distance_file, delimiter=',')[:,1:]
D = {label.replace("\"", "") : 
     {label.replace("\"", "") : element for label, element in zip(labels, row)} 
     for label, row in zip(labels, raw_data)}

main(graph_file, gml_filename, D=D)


# In[ ]:

import os
import networkx as nx
import numpy as np
from ghsom import main_no_labels as ghsom_main
import pickle
import shutil

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)
    
root_dir = "/home/david/Documents/ghsom"

data = "yeast_union_rel"
init = 1

for p in np.arange(0.1, 1, 0.1)[::-1]:
    
    print "p={}".format(p)
    
    os.chdir(root_dir)
    
    #ghsom parameters
    params = {'w': 0.0001,
             'eta': 0.001,
             'sigma': 1,
              'e_sg': p,
             'e_en': 10}
    
    map_file = '{}_communities_{}_{}'.format(data, p, init)
    
    if not os.path.isfile("{}.pkl".format(map_file)):
    
        #run ghsom and save output
        print "running GHSOM and saving to {}.pkl".format(map_file)
        G, map = ghsom_main(params, 'embedded_{}.gml'.format(data), init=init, lam=1000)
        print '\nnumber of communities detected: {}, saved map to {}'.format(len(map), map_file)
        save_obj((G, map), map_file)
    
    else:
        
        print "{}.pkl already exists, loading map".format(map_file)    
        #load output
        G, map = load_obj(map_file)

    #save results to file
    dir_name = "{}_communities_{}_{}".format(data, p, init)
    if not os.path.isdir(dir_name):
#         shutil.rmtree(dir_name)
#         print "deleted directory {}".format(dir_name)
    
        os.mkdir(dir_name)
        print 'made directory {}'.format(dir_name)

    os.chdir(dir_name)
    print "moved to {}".format(dir_name)
    
    #all genes
    all_genes_file = "all_genes.txt"
    with open(all_genes_file, 'w') as f:
        for n in G.nodes():
            f.write("{}\n".format(n))
    print "written {}".format(all_genes_file)
    
    #save shortest path matrix
    shortest_path = nx.floyd_warshall_numpy(map).astype(np.int)
    np.savetxt("shortest_path.csv", shortest_path, fmt='%i', delimiter=",")
    print 'written shortest path matrix'
    
    #save communities to file
    c = 0
    for n, d in map.nodes(data=True):
        ls = d['ls']
        with open('community_{}.txt'.format(c),'w') as f:
            for l in ls:
                f.write('{}\n'.format(l))
        print 'written community_{}.txt'.format(c)
        c += 1
    print


# In[ ]:



