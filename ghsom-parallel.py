
# coding: utf-8

# In[ ]:

from __future__ import division
import sys

import numpy as np
import networkx as nx
import sklearn.metrics as met

import matplotlib.pyplot as plt

from itertools import repeat

from Queue import Queue
from threading import Thread
from threading import current_thread

NUM_THREADS = 5
MIN_EXPANSION_SIZE = 10

##########################################################################################################################

#function to generate real valued som for graph input
def initialise_network(X):
    
    #network will be a one dimensional list
    network = nx.Graph()
    
    #initialise a network with just one neuron
    network.add_node(0)
    
    #assign a random vector in X to be the weight
    network.node[i]["v"] = X[np.random.randint(len(X))]
    
    #return network
    return network

##########################################################################################################################

# function to train SOM on given graph
def train_network(X, network, num_epochs, eta_0, sigma_0):
    
    #initial learning rate
    eta = eta_0
    
    #initial neighbourhood size
    sigma = sigma_0
    
    #list if all patterns to visit
    training_patterns = range(len(X))
    
    for e in range(num_epochs):
        
        #shuffle nodes
        np.random.shuffle(training_patterns)
        
        # iterate through N nodes of graph
        for i in trainingPatterns:
            
            #data point to consider
            x = X[i]
            
            #determine winning neuron
            closest_neuron = winning_neuron(x, network)
            
            # update weights
            update_weights(x, network, closest_neuron, eta, sigma)
            
        # drop neighbourhood
        sigma = sigma_0 * np.exp(-2 * sigma_0 * e / num_epochs);
        
##########################################################################################################################

# winning neuron
def winning_neuron(x, network):
    
    distances = {n: np.linalg.norm(x - d["v"]) for n, d in network.nodes(data=True)}
    
    return min(distances, key=distances.get)

##########################################################################################################################

# function to update weights
def update_weights(x, network, winning_neuron, eta, sigma):
        
    shortest_path_length = np.array([v for k, v in nx.shortest_path_length(network, source=winning_neuron).items()])
    
    weights = np.array([v for n, v in nx.get_node_attributes(network, "v").items()])
    weights += eta * np.exp(- shortest_path_length ** 2 / (2 * sigma ** 2)) * (x - weights)
    
    nx.set_node_attributes(network, "v", {n: v for n, v in zip(network.nodes(), weights)})

###########################################################################################################################
    
# assign nodes into clusters
def assign_nodes(names, X, network):
    
    distances = np.array([[np.linalg.norm(d["v"] - x) for n, d in network.nodes(data=True)] for x in X])
    min_distances = np.min(distances, axis=1)
    assignments = np.argmin(distances, axis=1)
    
    errors = np.array([np.mean(distances[assignments == n]) for n in network.nodes()])
    
    ##TODO delete nodes with error == NAN
    empty_neurons = np.where([np.isnan(errors)])[0]
    
    if empty_neurons.size > 0:
    
        ################################################DELETION####################################################

        #neighbours of deleted neurons
        neighbour_lists = np.array([network.neighbors(n) for n in empty_neurons])

        #remove the nodes
        network.remove_node_from(empty_neurons)
        
        ##connect separated components
        for neighbour_list in neighbour_lists:
            connect_components(network, neighbour_list)

        #########################################################################################################
    
    MQE = np.mean(errors)
    
    errors = {n: e for n, e in zip(network.nodes(), errors)}
    assignments = {n: [names[i] for i in np.where(assignments==n)] for n in network.nodes()}
    
    nx.set_node_attributes(network, "e", errors)
    nx.set_node_attributes(network, "ls", assignments)
    
    return MQE, len(empty_neurons)

##########################################################################################################################

def connect_components(network, neighbour_list):
    
    sub_network = network.subgraph(neighbour_list)
    
    connected_components = [c for c in nx.connected_components(sub_network)]
    number_of_connected_components = len(connected_components)
    
    for i in range(number_of_connected_components):
        
        for j in range(i + 1, number_of_connected_components):
            
            distances = np.array([[np.linalg.norm(n1["v"] - n2["v"]) for n2 in connected_components[j]] 
                                  for n1 in connected_components[i]])
            
            min_n1, min_n2 = np.unravel_index(distances.argmin(), distances.shape)
            
            network.add_edge(connected_components[i][min_n1], 
                            connected_components[j][min_n2])

##########################################################################################################################
            
##function to identify neuron with greatest error
def identify_error_unit(network):
    
    errors = nx.get_node_attributes(network, "e")
    
    return max(errors, key=errors.get)

##########################################################################################################################

def expand_network(named_X, network, error_unit):
    
    #id of new node
    id = max(network) + 1
    
    network.add_node(id)
    
    #v goes to random vector in range of error unit
    ls = network.node[error_unit]['ls']
    
    r = np.random.randint(len(ls))
    
    v = named_X[ls[r]]
    network.node[id]['v'] = v
        
    #identify neighbour pointing closet
    error_unit_neighbours = network.neighbors(error_unit)
    
    if len(error_unit_neighbours) == 0:
        
        ##add edge
        network.add_edge(error_unit, id)
        
    else:
        
        error_unit_vector = network.node[error_unit]["v"]
        
        ##find closest neighbour
        distances = {n: np.linalg.norm(v - n["v"]) for n in error_unit_neighbours}
        closest_neighbour = min(distances, key=distances.get)
        
        #connect to error unit and closest neighbour
        network.add_edge(closest_neighbour, id)
        network.add_edge(error_unit, id)
        
##########################################################################################################################
##########################################################################################################################

##GHSOM algorithm
def ghsom(ID, named_X, lam, eta, sigma, e_0, e_sg, e_en, q):
    
    #separate names and matrix of node embedding
    names, X = zip(*named_X.items())
    X = np.array(X)
    
    #create som for this neuron
    network = initialise_network(X)
    
    #train for lamda epochs
    train_network(X, network, lam, eta, sigma)

    #classify nodes and compute error
    MQE, num_deleted_neurons = assign_nodes(names, X, network)
    
    ##som growth phase
    #repeat until error is low enough
    while MQE > e_sg * e_0 and num_deleted_neurons < 3:
    
        #find neuron with greatest error
        error_unit = identify_error_unit(network)
        
        #expand network
        expand_network(named_X, network, error_unit)
        
        #train for lam epochs
        train_network(X, network, lam, eta, sigma)

        #calculate mean network error
        MQE, deleted_neurons = assign_nodes(names, X, network)
        num_deleted_neurons += deleted_neurons
        
        
    ##neuron expansion phase
    #iterate thorugh all neruons and find neurons with error great enough to expand
    for i, d in network.node(data=True):
        
        #unpack
        ls = d["ls"]
        e = d["e"]
        
        #check error
        if (e > e_en * e_0 and len(ls) > MIN_EXPANSION_SIZE and num_deleted_neurons < 3):
            
            id = "{}-{}".format(ID, str(i).zfill(2)) 
                
            sub_X = {k: named_X[k] for k in ls}
            
            print "submitted job: id={}, e={}".format(id, e)
            
            #add these parameters to the queue
            q.put((id, sub_X, lam, eta, sigma, e, e_sg, e_en))
    
    #return network
    return network, MQE

##########################################################################################################################
##########################################################################################################################

## get embedding TERRIBLE but staying
def get_embedding(G):
    
    #get number of niodes in the graph
    num_nodes = nx.number_of_nodes(G)
    
    #dimension of embedding
    d = 0
    while 'embedding'+str(d) in G.node[G.nodes()[0]]:
        d += 1
    
    #initialise embedding
    X = np.zeros((num_nodes, d))
    
    for i in range(num_nodes):
        for j in range(d):
            X[i,j] = G.node[G.nodes()[i]]['embedding'+str(j)]
    
    return X

####################################################################################################################

def worker(q, networks):
    
    while True:
        
        #unpack first element of queue
        #conatains all the para,eters for GHSOM
        ID, X, lam, eta, sigma, e_0, e_sg, e_en = q.get()
        
        #run GHSOM and return a network and MQE
        n, e = ghsom(ID, X, lam, eta, sigma, e_0, e_sg, e_en, q)
        
        #append result to networks list
        networks.append((ID, n, e))
        
        #mark task as done
        q.task_done()

def main(params, gml_filename, lam=10000):
    
    #network
    G = nx.read_gml(gml_filename)
    
    #embedding matrix
    X = get_embedding(G)
    
    #zip with names
    named_X = {k: v for k, v in zip(G.nodes(), X)}
    
    ##list of returned networks
    networks = []
    
    #initilise worker queue
    q = Queue()
    
    #initialise threads
    for i in range(NUM_THREADS):
        
        t = Thread(target=worker, args=(q, networks))
        t.setDaemon(True)
        t.start()
        
    ##initial MQE is variance of dataset
    m = np.mean(X, axis=0)
    MQE_0 = np.mean([np.linalg.norm(x - m) for x in X])
        
    #add initial layer of ghsom to queue
    q.put(("01", X, lam, params["eta"], params["sigma"], MQE_0, params["e_sg"], params["e_en"]))
    
    #finally wait until queue is empty and all tasks are done
    q.join()
    
    print "DONE"
    
    return G, networks


# In[ ]:

params = {'eta': 0.0001,
         'sigma': 1,
          'e_sg': 0.8,
         'e_en': 0.1}


# In[ ]:

# %%time 
get_ipython().magic(u'prun G, networks = main(params=params, gml_filename="embedded_yeast_uetz.gml", lam=100)')


# In[ ]:



