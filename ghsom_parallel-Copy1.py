
# coding: utf-8

# In[1]:

from __future__ import division

from sys import stdout

import numpy as np
import networkx as nx
import sklearn.metrics as met
from sklearn.metrics.pairwise import euclidean_distances

import matplotlib.pyplot as plt

from itertools import repeat

from Queue import Queue
from threading import Thread
from threading import current_thread

MIN_EXPANSION_SIZE = 10
MAX_DELETED_NEURONS = 3

#########################################################################################################################

##function to visualise graph
def visualise_graph(G, colours, layer):
        
    ## create new figure for graph plot
    fig, ax = plt.subplots()
    
    # graph layout
    pos = nx.spring_layout(G)
    
    #attributes in this graph
    attributes = np.unique([v for k, v in nx.get_node_attributes(G, "assigned_community_layer_{}".format(layer)).items()])

    # draw nodes -- colouring by cluster
    for i in range(min(len(colours), len(attributes))):
       
        node_list = [n for n in G.nodes() if G.node[n]["assigned_community_layer_{}".format(layer)] == attributes[i]]
        colour = [colours[i] for n in range(len(node_list))]
        
        nx.draw_networkx_nodes(G, pos, nodelist=node_list, node_color=colour)
        
    #draw edges
    nx.draw_networkx_edges(G, pos)

    # draw labels
    nx.draw_networkx_labels(G, pos)
    
    #title of plot
    plt.title('Nodes coloured by cluster, layer: {}'.format(layer))

    #show plot
    plt.show()

## visualise graph based on network clusters
def visualise_network(network, colours, layer):
    
    #num neurons in lattice
    num_neurons = len(network)

    ##create new figure for lattice plot
    fig, ax = plt.subplots()
    
    # graph layout
    pos = nx.spring_layout(network)

    # draw nodes -- colouring by cluster
    for i in range(len(colours)):
        nx.draw_networkx_nodes(network, pos, nodelist = [network.nodes()[i]], node_color = colours[i])

    #draw edges
    nx.draw_networkx_edges(network, pos)

    # draw labels
    nx.draw_networkx_labels(network, pos)
    
    #label axes
    plt.title('Neurons in lattice, layer: '+str(layer))
    
    #show lattice plot
    plt.show()

##########################################################################################################################

#function to generate real valued som for graph input
#three initial nodes
def initialise_network(ID, X, starting_nodes=3):
    
    #network will be a one dimensional list
    network = nx.Graph(ID = ID)
    
    #initialise a network with just one neuron
    network.add_nodes_from(range(1, starting_nodes + 1))
    
    #id of nodes
    for n in network.nodes():
        network.node[n]["ID"] = "{}-{}".format(ID, str(n).zfill(2))
        
    #connect nodes     
    for i in range(1, starting_nodes + 1):
        for j in range(i + 1, starting_nodes + 1):
            network.add_edge(i, j)
    
    #assign a random vector in X to be the weight
    V = X[np.random.randint(len(X), size=starting_nodes)]
    
#     print V
    
    #return network
    return network, V

#########################################################################################################################

def precompute_sigmas(sigma, num_epochs):
    
    return np.array([sigma * np.exp(-2 * sigma * e / num_epochs)
                     for e in range(num_epochs)])

##########################################################################################################################
##TODO
# function to train SOM on given graph
def train_network(X, network, V, num_epochs, eta_0, precomputed_sigmas):
    
    #initial learning rate
    eta = eta_0
    
    #list if all patterns to visit
    training_patterns = range(len(X))
    
    #shortest path matrix
    shortest_path = np.array(nx.floyd_warshall_numpy(network))
#     print "TRAINING"
#     print V
    
    for e in range(num_epochs):
        
        #shuffle nodes
        np.random.shuffle(training_patterns)
        
        sigma = precomputed_sigmas[e]
        
        # iterate through N nodes of graph
        for i in training_patterns:
            
            #data point to consider
            x = X[i]
            
            #determine winning neuron
            closest_neuron = winning_neuron(x, V)
            
            # update weights
#             V = update_weights(x, V, closest_neuron, shortest_path[closest_neuron], eta, pre_computed_sigmas[e])
            
            #weight update (vectorised)
            V += np.dot(np.diag(eta * np.exp(- shortest_path[closest_neuron] ** 2 / (2 * sigma * sigma))), (x - V))
            
#         stdout.write("\rTraining epoch: {}/{}".format(e, num_epochs))
            
#     print V
#     print
    return V
        
##########################################################################################################################

# winning neuron
def winning_neuron(x, V):
    
    distances = np.linalg.norm(x - V, axis=1)
    
    return distances.argmin()

##########################################################################################################################

# function to update weights
def update_weights(x, V, winning_neuron, shortest_path_length, eta, sigma):
    
    #weight update (vectorised)
    V += np.dot(np.diag(eta * np.exp(- shortest_path_length ** 2 / (2 * sigma ** 2))), 
                      (x - V))
    
    return V

########################################################################################################################   

# assign nodes into clusters
def assign_nodes(names, X, network, V):
    
    #distance from each datapoint (row) to each weight vector (column)
    distances = euclidean_distances(X, V)
    
    #minium distance for each datapoint
    min_distances = np.min(distances, axis=1)
    
    #index of column giving minimum distance
    arg_min_distances = np.argmin(distances, axis=1)
    
#     print "ARG MIN DISTANCES"
#     print arg_min_distances
    
    #nodes corresponding to minimum index (of length len(X))
    minimum_nodes = np.array([network.nodes()[n] for n in arg_min_distances])
    
#     print "MINIMUM NODES"
#     print minimum_nodes
    
    #list of neurons with no assignments
    empty_neurons = np.array([n for n in network.nodes() if n not in minimum_nodes])
    
    if empty_neurons.size > 0:
    
        ################################################DELETION####################################################

        #neighbours of deleted neurons
        neighbour_lists = np.array([network.neighbors(n) for n in empty_neurons])
        
#         print "DELETING NODES: {}".format(empty_neurons)
        
        #remove the nodes
        network.remove_nodes_from(empty_neurons)
        
        ##remove from V
        V = np.array([V[i] for i in range(len(V)) if i in arg_min_distances])
        
        #compute distances between all neurons in input space
        computed_neuron_distances = compute_euclidean_distances(network, V)
        
        ##connect separated components
        for neighbour_list in neighbour_lists:
            connect_components(network, neighbour_list, computed_neuron_distances)

        ############################################################################################################

    #array of errors
    errors = np.array([np.mean(min_distances[minimum_nodes == n]) for n in network.nodes()])
    
#     print "ERRORS"
#     print errors
    
    #compute MQE
    MQE = np.mean(errors)
    
    print "MQE={}, size of map={}".format(MQE, len(network))
    
    ##array of assignments
    assignments = np.array([np.array([names[i] for i in np.where(minimum_nodes == n)[0]]) for n in network.nodes()])
    
#     print "NAMES"
#     print names
    
#     print "ASSIGNMENTS"
#     print assignments
    
    #zip zith nodes
    errors = {n: e for n, e in zip(network.nodes(), errors)}
    assignments = {n: a for n, a in zip(network.nodes(), assignments)}
    
    nx.set_node_attributes(network, "e", errors)
    nx.set_node_attributes(network, "ls", assignments)
    
    return MQE, empty_neurons.size, V

##########################################################################################################################

def compute_euclidean_distances(network, V):
    
    distances = euclidean_distances(V)
    
    return {network.nodes()[i] : {network.nodes()[j] : distances[i, j] for j in range(len(distances[i]))}
           for i in range(len(distances))}
    
######################################################################################################################### 


def connect_components(network, neighbour_list, computed_neuron_distances):
    
    sub_network = network.subgraph(neighbour_list)
    
    connected_components = [sub_network.subgraph(c) for c in nx.connected_components(sub_network)]
    number_of_connected_components = len(connected_components)
    
    for i in range(number_of_connected_components):
        
        connected_component_1 = connected_components[i].nodes()
        
        for j in range(i + 1, number_of_connected_components):
            
            connected_component_2 = connected_components[j].nodes()
            
            distances = np.array([[computed_neuron_distances[n1][n2] for n2 in connected_component_2]
                                 for n1 in connected_component_1])
            
            min_n1, min_n2 = np.unravel_index(distances.argmin(), distances.shape)
            
            network.add_edge(connected_component_1[min_n1], 
                            connected_component_2[min_n2])

##########################################################################################################################
            
##function to identify neuron with greatest error
def identify_error_unit(network):
    
    errors = nx.get_node_attributes(network, "e")
    
    return max(errors, key=errors.get)

##########################################################################################################################

def expand_network(ID, named_X, network, V, error_unit):
    
    #v goes to random vector in range of error unit
    ls = network.node[error_unit]["ls"]    
    r = np.random.randint(len(ls))
    v = named_X[ls[r]]
    
    #zip nodes and distances
    distances = zip(network.nodes(), np.linalg.norm(V - v, axis=1))
        
    #identify neighbour pointing closet
    error_unit_neighbours = network.neighbors(error_unit)
    
    
    #id of new node
    new_node = max(network) + 1
    
    #add new node to map
    network.add_node(new_node)
    
    ##id
    network.node[new_node]["ID"] = "{}-{}".format(ID, str(new_node).zfill(2))
    
    #add edges to map
    
    #connect error unit and new node
    network.add_edge(error_unit, new_node)
    
    if len(error_unit_neighbours) > 0:
        
        ##find closest neighbour
        distances = {n: v for n, v in distances if n in error_unit_neighbours}
        closest_neighbour = min(distances, key=distances.get)
        
        #connect to error unit and closest neighbour
        network.add_edge(closest_neighbour, new_node)
        
    #add v to V
    V = np.vstack([V, v])   
    
#     print V
    
    return V
        
##########################################################################################################################
##########################################################################################################################

##GHSOM algorithm
def ghsom(ID, named_X, lam, eta, sigma, e_0, e_sg, e_en, q):
    
    print "MQE_0={}, growth target={}".format(e_0, e_0 * e_sg)
    
    #separate names and matrix of node embedding
    names, X = zip(*named_X.items())
    names = np.array(names)
    X = np.array(X)
    
    #create som for this neuron
    network, V = initialise_network(ID, X)
    
    #precompute sigmas
    precomputed_sigmas = precompute_sigmas(sigma, lam)
    
    #train for lamda epochs
    V = train_network(X, network, V, lam, eta, precomputed_sigmas)
    
    #classify nodes and compute error
    MQE, num_deleted_neurons, V = assign_nodes(names, X, network, V)
    
    ##som growth phase
    #repeat until error is low enough
    while MQE > e_sg * e_0 and num_deleted_neurons < MAX_DELETED_NEURONS:
        
        #find neuron with greatest error
        error_unit = identify_error_unit(network)
        
        #expand network
        V = expand_network(ID, named_X, network, V, error_unit)
        
        #train for lam epochs
        V = train_network(X, network, V, lam, eta, precomputed_sigmas)

        #calculate mean network error
        MQE, deleted_neurons, V = assign_nodes(names, X, network, V)
        num_deleted_neurons += deleted_neurons
        
    print "growth terminated, MQE: {}, target: {}, number of deleted neurons: {}".format(MQE, e_0 * e_sg, num_deleted_neurons)
        
    ##neuron expansion phase
    #iterate thorugh all neruons and find neurons with error great enough to expand
    for _, d in network.nodes(data=True):
        
        #unpack
        node_id = d["ID"]
        ls = d["ls"]
        e = d["e"]
        
        #check error
        if (e > e_en * e_0 and len(ls) > MIN_EXPANSION_SIZE and num_deleted_neurons < MAX_DELETED_NEURONS):
            
#             id = "{}-{}".format(ID, node_id) 
                
            sub_X = {k: named_X[k] for k in ls}
            
            print "submitted job: ID={}, e={}, number of nodes={}".format(node_id, e, len(ls))
            
            #add these parameters to the queue
            q.put((node_id, sub_X, lam, eta, sigma, e, e_sg, e_en))
    
    #return network
    return network, MQE

##########################################################################################################################
##########################################################################################################################

def label_nodes(G, networks):
    
    for _, network, _ in networks: 
        
        for _, d in network.nodes(data=True):
            
            community = d["ID"]
            layer = community.count("-")
            assignment_string = "assigned_community_layer_{}".format(layer)
            
            for node in d["ls"]:
                
                G.node[node][assignment_string] = community


##########################################################################################################################

def NMI_one_layer(G, label, layer):
    
    #actual community for this layer
    actual_community_labels = np.array([v for k, v in nx.get_node_attributes(G, label).items()])
    
    #predicted communitiy for this layer
    predicted_community_labels = np.array([v for k, v in nx.get_node_attributes(G, 
         "assigned_community_layer_{}".format(layer)).items()])

    print actual_community_labels
    print predicted_community_labels
    
    return met.normalized_mutual_info_score(actual_community_labels, predicted_community_labels)

def NMI_all_layers(G, labels):
    
    return np.array([NMI_one_layer(G, labels[i], i + 1) for i in range(len(labels))])

##########################################################################################################################

## get embedding TERRIBLE but staying
def get_embedding(G):
    
#     #get number of niodes in the graph
#     num_nodes = nx.number_of_nodes(G)
    
#     #dimension of embedding
#     dim = 0
#     while 'embedding'+str(dim) in G.node[G.nodes()[0]]:
#         dim += 1
    
#     #initialise embedding
#     X = np.array([[d["embedding{}".format(j)] for j in range(dim)] for n, d in G.nodes(data=True)])

    return np.array([v for k, v in nx.get_node_attributes(G, "embedding").items()])


##########################################################################################################################

def process_job(q, networks):
    
    #unpack first element of queue
    #contains all the para,eters for GHSOM
    ID, X, lam, eta, sigma, e_0, e_sg, e_en = q.get()

    #run GHSOM and return a network and MQE
    n, e = ghsom(ID, X, lam, eta, sigma, e_0, e_sg, e_en, q)

    #append result to networks list
    networks.append((ID, n, e))

    #mark task as done
    q.task_done()

def worker(q, networks):
    
    #continually poll queue for jobs 
    while True:
        process_job(q, networks)

def main(params, filename, lam=10000, num_threads=1):
    
    #network
    G = nx.read_gpickle(filename)
    
    #embedding matrix
    X = get_embedding(G)
    
    #zip with names
    named_X = {k: v for k, v in zip(G.nodes(), X)}
    
    ##list of returned networks
    networks = []
    
    #initilise worker queue
    q = Queue()
    
    ##initial MQE is variance of dataset
    m = np.mean(X, axis=0)
    MQE_0 = np.mean([np.linalg.norm(x - m) for x in X])
    
    #add initial layer of ghsom to queue
    q.put(("01", named_X, lam, params["eta"], params["sigma"], MQE_0, params["e_sg"], params["e_en"]))
    
    if num_threads > 1:
    
        #initialise threads
        for i in range(num_threads):

            t = Thread(target=worker, args=(q, networks))
            t.setDaemon(True)
            t.start()

        #finally wait until queue is empty and all tasks are done
        q.join()
        
    else :
        
        #single thread
        while not q.empty():
            process_job(q, networks)
    
    print "DONE"
    
    return G, networks


# In[68]:

# params = {'eta': 0.001,
#          'sigma': 1,
#           'e_sg': 0.5,
#          'e_en': 1.0}


# In[72]:

# G, networks = main(params=params, filename="embedded_benchmark.gpickle", num_threads=5, lam=1000)


# In[70]:

# label_nodes(G, networks)


# In[73]:

# NMI_all_layers(G, labels=["firstlevelcommunity"])


# In[74]:

# _, network, _ = networks[0]
# colours = np.random.rand(len(network), 3)


# In[75]:

# visualise_graph(G=G, colours=colours, layer=1)

