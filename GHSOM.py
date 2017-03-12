
# coding: utf-8

# In[1]:

import numpy as np
import networkx as nx
import som_functions as som
import math
import sklearn.metrics as met
import sklearn.manifold as man
from time import time


# In[48]:

##GHSOM algorithm

#G: graph
#lam: lambda -- the number of epochs to train before assessing error
#eta: learning rate
#sigma: initial neighbourhood range
#e_0: error of previous layer
#e_sg: error must reduce by this much for growth to stop
#e_en: error must be greater than this to expand neuron
#layer: layer of som
#n: desired initial number of neurons
#m: desired number of neurons in the next layer
def ghsom(G, lam, w, eta, sigma, e_0, e_sg, e_en, layer, n, m):
    
    #embedding
    X = dsd_embedding(G)
    
    #save embedding to graph
    set_embedding(G, X)
    
    print('embedded graph')
    
    #number of nodes in G
    num_nodes = nx.number_of_nodes(G)
    
    ##number of training patterns to visit
    N = min(num_nodes, 100)
#     N = num_nodes
    
    #create som for this neuron
    network = som.initialise_network(X, n, w)
    
    ##inital training phase
    
    #train for lam epochs
    network = som.train_network(X, network, lam, eta, sigma, N)

    #classify nodes
    network = som.assign_nodes(G, X, network, layer)

    #calculate mean network error
    network, MQE = som.update_errors(network)
    
    ##som growth phase
    #repeat until error is low enough
    while MQE > e_sg * e_0:
#     for i in range(1):
#     if layer > 0:
    
        #find neuron with greatest error
        error_unit = som.identify_error_unit(network)
        
        #expand network
        network = som.expand_network(network, error_unit)
        
        #train for l epochs
        network = som.train_network(X, network, lam, eta, sigma, N)

        #classify nodes
        network = som.assign_nodes(G, X, network, layer)

        #calculate mean network error
        network, MQE = som.update_errors(network)
        
        print('ghsom has expanded som',layer,'error',MQE)
        
    print('ghsom has terminated expansion',layer)
    print('error',MQE)
    
    #recalculate error after neuron expansion
    MQE = 0
    
    ##neuron expansion phase
    #iterate thorugh all neruons and find neurons with error great enough to expand
    for i in network.nodes():
        
        #unpack
        ls = network.node[i]['ls']
        e = network.node[i]['e']
        
        #check error
        if e > e_en * e_0 or e_0 == math.inf:
#         if layer < len(labels) and len(ls) > 0:

            if e_0 == math.inf:
                e_0 = e
        
            #subgraph
            H = G.subgraph(ls)
            
            #recursively run algorithm to create new network for subgraph of this neurons nodes
            n, e = ghsom(H, lam, w, eta, sigma, e_0, e_sg, e_en, layer + 1, m, m)
            
            #repack
            network.node[i]['e'] = e
            network.node[i]['n'] = n
            
            print('ghsom has built new layer',layer+1)
            
        #increase overall network error
        MQE += e
    
    #mean MQE
    MQE /= nx.number_of_nodes(network)
    
    #return network
    return network, MQE


# In[3]:

##visualise network
def visualise(G, network, neurons_in_each_layer, layer):
    
    ##randomly generate colours for plotting
    colours = np.random.rand(len(network), 3)
    
    ##visualise in ghsom function
    som.visualise_graph(G, colours, neurons_in_each_layer, layer)
    som.visualise_network(network, colours, layer)
    
    for i in network.nodes():
        
        #unpack
        l = network.node[i]['ls']
        n = network.node[i]['n']
        
        if len(n) > 0:

            H = G.subgraph(l)

            visualise(H, n, neurons_in_each_layer, layer + 1)


# In[4]:

##function to recursively label nodes in graph
def label_graph(G, network, layer, neuron_count):
    
    #number of neurons in network
    num_neurons = len(network)
    
    for i in network.nodes():
        
        #unpack
        l = network.node[i]['ls']
        
        for node in l:
            G.node[node]['community'+str(layer)] = neuron_count[layer]
            
        n = network.node[i]['n']
            
        if len(n) > 0: 
            
            H = G.subgraph(l)
            
            G = label_graph(G, n, layer + 1, neuron_count)
            
        neuron_count[layer] += 1
            
    return G


# In[5]:

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


# In[6]:

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


# In[7]:

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


# In[8]:

##function to calculate community detection error given a generated benchmark graph
def mutual_information(G, labels):
    
    #number of layers of cluser
    num_layers = len(labels)
    
    #initialise scores
    scores = [0 for l in range(num_layers)]
    
    #iterate over all levels of labels
    for i in range(num_layers):
    
        #assigned first layer community
        actual_community = nx.get_node_attributes(G, labels[num_layers - i - 1])

        #predicted first layer community
        predicted_community = nx.get_node_attributes(G, 'community'+str(i + 1))

        #labels for first layer of community
        labels_true = [v for k,v in actual_community.items()]
        labels_pred = [v for k,v in predicted_community.items()]
        
        if len(labels_pred) == 0:
            continue

        #mutual information to score classifcation
        scores[i] = met.normalized_mutual_info_score(labels_true, labels_pred)
#         scores[i] = met.adjusted_mutual_info_score(labels_true, labels_pred)
    
    #return
    return scores


# In[9]:

##function to compute G = (I - P + W)^-1
#N: normalised laplacian
#d: list of degrees
#v: eigenvector corresponding to smallest eigenvalue
def compute_gr(N, B, d, v):
    
    #number of nodes in the graph
    n = len(N)
    
    #initialise Gr
    Gr = np.zeros((n, n))
    
    for i in range(n):
        
        #b2
        b2 = (np.transpose(v) @ B[:,i]) * v
        
        #b1 
        b1 = B[:,i] - b2
        
        #x1
        x1 = np.linalg.lstsq(N, b1)[0]
        
        #add to Gr
        Gr[:,i] = np.diag(d ** 0.5) @ (x1 + b2)
    
    #return Gr
    return Gr


# In[10]:

##function to compute DSD matrix using AMG method
#N: laplacian matrix
#d: list of degrees
#v: smallest eigenvector
def dsd(N, d, v):
    
    #number of nodes in G
    n = len(N)
    
    #initialize dsd matrix
    dsd = np.zeros((n, n))
    
    #B
    B = np.diag(d ** -0.5) @ np.identity(n)
    
    #compute G
    G = compute_gr(N, B, d, v)
    
    print('computed greens matrix')
    
    #compute dsd for each pair
    for i in range(n):
        for j in range(i, n):
            
            #compute distance
            dis = np.linalg.norm(np.transpose(B[:,i] - B[:,j]) @ G, ord=1)
            
            #add to dsd matrix
            dsd[i,j] = dis
            dsd[j,i] = dis
    
    #return
    return dsd


# In[11]:

##function to return all connected components
def connected_components(G):
    
    #number of nodes in the graph
    n = nx.number_of_nodes(G)
    
    #adjacency matrix
    A = nx.adjacency_matrix(G).toarray()
    
    #list of degrees
    deg = np.sum(A, axis=1)
    
    #normalised graph laplacian
    N = np.identity(n) - np.diag(deg ** -0.5) @ A @ np.diag(deg ** -0.5)
    
    #eigen decompose normalised laplacian
    l, U = np.linalg.eigh(N)
    
    ##sort eigenvalues (and eigenvectors) into ascending order
    idx = l.argsort()
    l = l[idx]
    U = U[:,idx]
    
    connected_components = len(l[l<1e-12])
    print('number of connected components',connected_components)
    
    connected_graphs = []
    
    for i in range(connected_components):
        
        ids = U[:,i].nonzero()[0]
        
        connected_graphs.append(G.subgraph([G.nodes()[n] for n in ids]))
        
    return connected_graphs


# In[41]:

##dsd embedding
def dsd_embedding(G):
    
    #number of nodes in the graph
    n = nx.number_of_nodes(G)
    
    #adjacency matrix
    A = nx.adjacency_matrix(G).toarray()
    
    #list of degrees
    deg = np.sum(A, axis=1)
    
    #normalised graph laplacian
    N = np.identity(n) - np.diag(deg ** -0.5) @ A @ np.diag(deg ** -0.5)
    
    print('constructed normalised laplacian')
    
    #eigen decompose normalised laplacian
    l, U = np.linalg.eigh(N)
    
    ##sort eigenvalues (and eigenvectors) into ascending order
    idx = l.argsort()
    l = l[idx]
    U = U[:,idx]
    
    #compute dsd matrix as metric
    D = dsd(N, deg, U[:,0]) 
    
    print('computed dsd matrix')
    
    #centreing matrix
    C = np.identity(n) - np.ones((n, n)) / n
    
    #similarity matrix
    K = - 1/2 * C @ D ** 2 @ C
    
    #eigen decompose K
    l, U = np.linalg.eigh(K)
    
    ##sort eigenvalues (and eigenvectors) into ascending order
    idx = l.argsort()[::-1]
    l = l[idx]
    U = U[:,idx]
    
    #sum of all eigen values
    s = sum(l)
    
    #estimate the number of dimensions to keep
    k = len(l)
    var = 1
    
    while var > 0.95:
        var = sum(l[:k]) / s
        k -= 1
    
    k += 1
   
    
    ##TODO
#     k = 3
    print('reduced dimension of data',k)
    
    #position matrix
    X = U[:,:k] @ np.diag(l[:k] ** 0.5)
    
    return X


# In[13]:

##save embedding to graph
def set_embedding(G, X):
    
    #get number of niodes in the graph
    num_nodes = nx.number_of_nodes(G)
    
    #iterate over a;; the nodes and save their embedding
    for i in range(num_nodes):
        G.node[G.nodes()[i]]['embedding'] = X[i]        


# In[14]:

##get embedding from graph
def get_embedding(G):
    #get the number of nodes in the graph
    num_nodes = nx.number_of_nodes(G)
    
    # get the embeddings
    embeddings = nx.get_node_attributes(G, 'embedding')
    
    #append all embeddings into one array
    X = np.array([v for k,v in embeddings.items()])
    
    #return the array
    return X


# In[37]:

# generate graph
# G = nx.karate_club_graph()
# G = nx.read_gml('football.gml')
# G = nx.read_gml('netscience.gml')
# G = nx.read_gml('dolphins.gml')
# G = nx.read_gml('karate.gml')
# G = benchmark_graph('bin_network.dat', 'community.dat')
# G = benchmark_hierarchical_graph('network.dat', 
#                                  'community_first_level.dat', 
#                                  'community_second_level.dat')
# G = benchmark_hierarchical_graph('literature_network.dat', 
#                                  'literature_community_first_level.dat', 
#                                  'literature_community_second_level.dat')
G = benchmark_hierarchical_graph('literature_network_32.dat', 
                                 'literature_community_first_level_32.dat', 
                                 'literature_community_second_level_32.dat')
# G = benchmark_hierarchical_graph('literature_network_double.dat', 
#                                  'literature_community_first_level_double.dat', 
#                                  'literature_community_second_level_double.dat')
print('loaded G')

# labels = ['club']
# labels = ['firstlevelcommunity']
labels = ['firstlevelcommunity','secondlevelcommunity']

#divide graph into connected components
connected_graphs = connected_components(G)
print('calculated number of connected components')


# In[ ]:

##start time
start = time()

#number of epochs to train == lambda
lam = 1000

#weights are intialised based on a uniform random distribution between +- this value
w = 1e-4

#initial learning rate
eta_0 = 1e-3

#initial neighbourhood size
sigma_0 = 1

#stop growth of current layer
e_sg = 0.8

#error must be greater than this times previous error for expansion
e_en = 0.9

#layer of GHSOM
layer = 0;

#desired number of initial neurons
n = 1

#desired number of initial neurons for all other layers
m = 2

for H in connected_graphs:

    #run ghsom algorithm
    network, MQE = ghsom(H, lam, w, eta_0, sigma_0, math.inf, e_sg, e_en, layer, n, m)

    print('ghsom algorithm has terminated')
    print('mean network error:',MQE)

    #label graph
    neurons = np.zeros(50, dtype=np.int)
    G = label_graph(G, network, layer, neurons)

    ##visualise
    neurons = np.zeros(50, dtype=np.int)
    visualise(G, network, neurons, 0)

    ##calculate error
    mi_score = mutual_information(G, labels)
    print('normalised mutual information score',mi_score)

##time taken
print('time taken',(time() - start))


# In[80]:

import matplotlib.pyplot as plt

##save network to gml file
N = nx.Graph();

N, num_nodes = build_connected_som(N, network, 0)

nx.draw_networkx(N)
plt.show()


# In[43]:

def build_connected_som(N, network, num):
    
    num_nodes = num
    
    for n in range(len(network.nodes())):

        num_nodes += 1
        
        node1 = network.nodes()[n]
        
        N.add_edge(num, num_nodes)
        
        for m in range(n):
            
            node2 = network.nodes()[m]
            
            if network.has_edge(node1, node2):
                
                N.add_edge(num_nodes, num + m + 1)     
        
        som = network.node[node1]['n']
        
        if len(som) > 0:
            
            N, num_nodes = build_connected_som(N, som, num_nodes)

        
    return N, num_nodes

