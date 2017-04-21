
# coding: utf-8

# In[1]:

from __future__ import division
import sys

import numpy as np
import networkx as nx
import sklearn.metrics as met

import matplotlib.pyplot as plt

from itertools import repeat

from Queue import Queue
from threading import Thread


NUM_THREADS = 5
MIN_EXPANSION_SIZE = 10

#function to generate real valued som for graph input
def initialise_network(X, num_neurons):
    
    #network will be a one dimensional list
    network = nx.Graph()
    
    #number of data points
    num_nodes = len(X) 
    
    #dimension of data in X
    d = len(X[0])
    
    #regular lattice
    lattice_size = np.floor(np.sqrt(num_neurons))
    
    for i in range(num_neurons):
        
        ##position
        network.add_node(i)
        
        ##weight    
        r = np.random.randint(num_nodes)
        network.node[i]['v'] = X[r]

        ##list of closest nodes
        network.node[i]['ls'] = []

        ##error of neuron
        network.node[i]['e'] = 0
        
        #connections
        if i % lattice_size > 0:
            #horizontal connection
            network.add_edge(i, i - 1)
        
        if i >= lattice_size:
            #vertical connection
            network.add_edge(i, i - lattice_size)
            
            if i % lattice_size < lattice_size - 1:
                #diagonal connection
                network.add_edge(i, i - lattice_size + 1)
                
    
    #return network
    return network


# function to train SOM on given graph
def train_network(X, network, num_epochs, eta_0, sigma_0, N, layer, MQE, target, num_deleted_neurons):
    
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
        for i in range(N):
            
            #data point to consider
            x = X[training_patterns[i]]
            
            #determine winning neuron
            closest_neuron = winning_neuron(x, network)
            
            # update weights
            update_weights(x, network, closest_neuron, eta, sigma)
            
        # drop neighbourhood
        sigma = sigma_0 * np.exp(-2 * sigma_0 * e / num_epochs);

# winning neuron
def winning_neuron(x, network):

    # minimum distance so far
    min_dist = np.inf
    winning_neuron = []
    
    # iterate through network
    for i in network.nodes():
            
        #unpack network
        v = network.node[i]['v']

        #distance between input vector and neuron weight
        distance = np.linalg.norm(x - v)

        # if we have a new closest neuron
        if distance < min_dist:
            min_dist = distance
            winning_neuron = i
    
    #return
    return winning_neuron

# function to update weights
def update_weights(x, network, winning_neuron, eta, sigma):
    
    # iterate through all neurons in network
    for i in network.nodes():
        
        #unpack
        v = network.node[i]['v']

        #new v -- move along shortest path by move distance
        v += eta * neighbourhood(network, i, winning_neuron, sigma) * (x - v)

        #save to network
        network.node[i]['v'] = v

# neighbourhood function
def neighbourhood(network, r, winning_neuron, sigma):
    
    return np.exp(-(nx.shortest_path_length(network, r, winning_neuron)) ** 2 / (2 * sigma ** 2))

# assign nodes into clusters
def assign_nodes(G, X, network, layer):
    
    #number of neurons in network
    num_neurons = nx.number_of_nodes(network)
    
    #number of nodes
    num_nodes = nx.number_of_nodes(G)
    
    # clear existing closest node list
    for i in network.nodes():
        
        network.node[i]['ls'] = []
        network.node[i]['e'] = 0
    
    # assign colour to each node
    for n in range(num_nodes):
        
        #data point to assign
        x = X[n]
    
        #intialise distance to be infinity
        min_distance = np.inf
        
        #closest reference vector to this ndoe
        closest_ref = []

        # find which neuron's referece vector this node is closest to
        for i in network.nodes():
                
            #unpack network
            v = network.node[i]['v']

            # calculate distance to that reference vector
            d = np.linalg.norm(x - v)

            if d < min_distance:
                min_distance = d
                closest_ref = i
        
        #unable to find closest ref
        if closest_ref == []:
            continue
        
        #add node to closest nodes list
        network.node[closest_ref]['ls'].append(G.nodes()[n])
        
        #increase e by distance
        network.node[closest_ref]['e'] += min_distance

##function to return lattice grid of errors
def update_errors(network, num_deleted_neurons):
    
    #mean network error
    mqe = 0;
    
    #neurons with assigned nodes
    num_neurons = 0
    
    #iterate over all neurons and average distance
    for i in network.nodes():
            
        #unpack network
        e = network.node[i]['e']
        ls = network.node[i]['ls']
        
        if len(ls) == 0:
            delete_node(network, i)
            num_deleted_neurons += 1
            continue
            
        num_neurons += 1

        #divide by len(ls) for mean
        e /= len(ls)

        #sum total errors
        mqe += e
        
        #save error to network
        network.node[i]['e'] = e        
    
    #mean
    mqe /= num_neurons
    
    return mqe, num_deleted_neurons

def connect_closest_neurons(network, s1, s2):
    
    min_dist = np.inf
    
    c1 = []
    c2 = []
    
    for i in s1:
        
        v1 = network.node[i]['v']
        
        for j in s2:
            
            v2 = network.node[j]['v']
            
            d = np.linalg.norm(v1 - v2)
            
            if d < min_dist:
                
                min_dist = d
                c1 = i
                c2 = j
            
    ##connect
    network.add_edge(c1, c2)
    
def intersection(a, b):
     return list(set(a) & set(b))

def delete_node(network, n):
    
    neighbours = network.neighbors(n)
    
    network.remove_node(n)
    
    if not neighbours:
        return
    
    components = [c for c in nx.connected_components(network)]
    
    for i in range(len(components)):
        
        conn_neigh_1 = intersection(components[i], neighbours)
        
        for j in range(i + 1, len(components)):
            
            conn_neigh_2 = intersection(components[j], neighbours)
            
            ##make connection between closest neurons in input space
            connect_closest_neurons(network, conn_neigh_1, conn_neigh_2)
            
    
            
##function to identify neuron with greatest error
def identify_error_unit(network):
    
    #initial value for maximum error found
    max_e = 0
    
    #initial index to return
    error_node = []
    
    for i in network.nodes():
        
        #unpack
        e = network.node[i]['e']
        
        #check if this unit has greater error than maximum 
        if e > max_e:
            
            max_e = e
            error_node = i
            
    #return id of unit with maximum error
    return error_node

def get_vector(node):
    
    d = 0
    
    while 'embedding'+str(d) in node:
        d += 1
    
    v = np.zeros(d)
    
    for i in range(d):
        v[i] = node['embedding{}'.format(i)]
        
    return v

def expand_network(G, network, error_unit):
    
    #id of new node
    id = max(network) + 1
    
    network.add_node(id)
    
    #v goes to random vector in range of error unit
    ls = network.node[error_unit]['ls']
    r = np.random.randint(len(ls))
    node = G.node[ls[r]]
    v = get_vector(node)
    network.node[id]['v'] = v
    
    ##list of closest nodes
    ls = []
    network.node[id]['ls'] = ls

    ##error of neuron
    e = 0
    network.node[id]['e'] = e

    ##som for neuron
    n = []
    network.node[id]['n'] = n
    
    #connections to other neurons
        
    #identify neighbour pointing furthest away
    error_unit_neighbours = network.neighbors(error_unit)
    
    if len(error_unit_neighbours) == 0:
        
        ##add edge
        network.add_edge(error_unit, id)
        
    else:
        
        ##find closest neighbour
        n = closest_neuron(network, error_unit, error_unit_neighbours)
        
        #connect to error unit and closest neighbour
        network.add_edge(n, id)
        network.add_edge(error_unit, id)

##function to find neuron pointing furthest away in list
def furthest_neuron(network, error_unit, ls):
    
    vi = network.node[error_unit]['v']
    max_dist = -np.inf
    
    #neighbour id
    furthest_node = []

    #iterate through neighbours
    for i in ls:

        #unpack neighbour
        v = network.node[i]['v']

        #distance in input space
        d = np.linalg.norm(v - vi)

        #is d > max_dist?
        if d > max_dist:

            max_dist = d
            furthest_node = i
            
    return furthest_node

def closest_neuron(network, error_unit, ls):
    
    vi = network.node[error_unit]['v']
    min_dist = np.inf
    
    #neighbour id
    closest_node = []

    #iterate through neighbours
    for i in ls:

        #unpack neighbour
        v = network.node[i]['v']

        #distance in input space
        d = np.linalg.norm(v - vi)

        #is d > max_dist?
        if d < min_dist:

            min_dist = d
            closest_node = i
            
    return closest_node

##GHSOM algorithm
def ghsom(ID, G, lam, eta, sigma, e_0, e_sg, e_en, init, layer, q):
    
    #embedding
    X = get_embedding(G)
    
    #number of nodes in G
    num_nodes = nx.number_of_nodes(G)
    
    ##number of training patterns to visit
    N = num_nodes
    
    if layer == 0:
        ini = 1
    else:
        ini = init
    
    #create som for this neuron
    network = initialise_network(X, ini)
    
    ##inital training phase
    
    #mean quantization error
    MQE = np.inf
    
    #number of neurons deleted from the map
    num_deleted_neurons = 0
    
    #train for lam epochs
    train_network(X, network, lam, eta, sigma, N, layer, MQE, e_sg * e_0, num_deleted_neurons)

    #classify nodes
    assign_nodes(G, X, network, layer)

    #calculate mean network error
    MQE, num_deleted_neurons = update_errors(network, num_deleted_neurons)
    
    ##som growth phase
    #repeat until error is low enough
    while MQE > e_sg * e_0 and num_deleted_neurons < 3:
    
        #find neuron with greatest error
        error_unit = identify_error_unit(network)
        
        #expand network
        expand_network(G, network, error_unit)
        
        #train for l epochs
        train_network(X, network, lam, eta, sigma, N, layer, MQE, e_sg * e_0, num_deleted_neurons)

        #classify nodes
        assign_nodes(G, X, network, layer)

        #calculate mean network error
        MQE, num_deleted_neurons = update_errors(network, num_deleted_neurons)
        
    ##neuron expansion phase
    #iterate thorugh all neruons and find neurons with error great enough to expand
    for i in range(len(network)):
        
        #unpack
        ls = network.node[network.nodes()[i]]['ls']
        e = network.node[network.nodes()[i]]['e']
        
        #check error
        if (e > e_en * e_0 and len(ls) > MIN_EXPANSION_SIZE and num_deleted_neurons < 3) or e_0 == np.inf:

            if e_0 == np.inf:
                e_0 = e
                
                id = "01"
            else:
                id = "{}-{}".format(ID, str(i).zfill(2)) 
        
            #subgraph
            H = G.subgraph(ls)
            
            #add these parameters to the queue
            print "added job (ID={}, e={}, layer={}) to queue".format(id, e, layer + 1)
            q.put((id, H, lam, eta, sigma, e, e_sg, e_en, init, layer + 1))
    
    #return network
    return network, MQE


def unassign_all_nodes(G, labels):

    #number of layers of communities
    num_layers = len(labels)

    for l in range(num_layers):

        nx.set_node_attributes(G, 'community'+str(l), 'unassigned')

##function to recursively label nodes in graph
def label_graph(G, network, layer, neuron_count):
    
    for i in network.nodes():
        
        #unpack
        l = network.node[i]['ls']
        
        for node in l:
            G.node[node]['community'+str(layer)] = neuron_count[layer]
            
        n = network.node[i]['n']
            
        if len(n) > 0: 
            
            H = G.subgraph(l)
            
            label_graph(H, n, layer + 1, neuron_count)
            
        neuron_count[layer] += 1


##function to calculate community detection error given a generated benchmark graph
def mutual_information(G, labels):
    
    #number of layers of communities
    num_layers = len(labels)
    
    #initialise scores
    scores = np.zeros(num_layers)
    
    #iterate over all levels of labels
    for i in range(num_layers):
    
        #assigned first layer community
        actual_community = nx.get_node_attributes(G, labels[num_layers - i - 1])

        #predicted first layer community
        predicted_community = nx.get_node_attributes(G, 'community'+str(i + 1))

        #only retrieve labels for assigned nodes
        labels_true = [v for k,v in actual_community.items() if k in predicted_community]
        labels_pred = [v for k,v in predicted_community.items() if k in actual_community]
        
        if len(labels_pred) == 0:
            continue
            
        #mutual information to score classifcation -- scale by number of assigned nodes out of all nodes
        score = met.normalized_mutual_info_score(labels_true, labels_pred) * len(labels_pred) / len(actual_community)
        scores[i] = score
    
    #return
    return scores 


## get embedding
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

# def modularity(G, H):
    
#     #number of links in communitiy H
#     l_s = nx.number_of_edges(H)
    
#     #total degree of communitiy H
#     d_s = np.sum(list(H.degree().values()))
    
#     #number of links in G
#     L = nx.number_of_edges(G)
        
#     #modularitiy
#     Q = l_s / L - (d_s / (2 * L)) ** 2
    
#     return Q

# ##function to visualise graph
# def visualise_graph(G, colours, layer):
        
#     ## create new figure for graph plot
#     fig, ax = plt.subplots()
    
#     # graph layout
#     pos = nx.spring_layout(G)
    
#     #attributes in this graph
#     attributes = np.unique([v for k,v in nx.get_node_attributes(G, 'community'+str(layer)).items()])

#     # draw nodes -- colouring by cluster
#     for i in range(min(len(colours), len(attributes))):
       
#         node_list = [n for n in G.nodes() if G.node[n]['community'+str(layer)] == attributes[i]]
#         colour = [colours[i] for n in range(len(node_list))]
        
#         nx.draw_networkx_nodes(G, pos, nodelist=node_list, node_color=colour)
        
#     #draw edges
#     nx.draw_networkx_edges(G, pos)

#     # draw labels
#     nx.draw_networkx_labels(G, pos, )
    
#     #title of plot
#     plt.title('Nodes coloured by cluster, layer: '+str(layer))

#     #show plot
#     plt.show()

# ## visualise graph based on network clusters
# def visualise_network(network, colours, layer):
    
#     #num neurons in lattice
#     num_neurons = len(network)

#     ##create new figure for lattice plot
#     fig, ax = plt.subplots()
    
#     # graph layout
#     pos = nx.spring_layout(network)

#     # draw nodes -- colouring by cluster
#     for i in range(len(colours)):
#         nx.draw_networkx_nodes(network, pos, nodelist = [network.nodes()[i]], node_color = colours[i])

#     #draw edges
#     nx.draw_networkx_edges(network, pos)

#     # draw labels
#     nx.draw_networkx_labels(network, pos)
    
#     #label axes
#     plt.title('Neurons in lattice, layer: '+str(layer))
    
#     #show lattice plot
#     plt.show()

# ###evaluate fitness
# def fitness(eta, sigma, e_sg, e_en, gml_filename, labels, init, lam):
    
#     G = nx.read_gml(gml_filename)
#     labels = labels.split(',')

#     #start layer
#     layer = 0
    
#     #run ghsom algorithm
#     network, MQE = ghsom(G, lam, eta, sigma, np.inf, e_sg, e_en, init, layer)

#     #label graph
#     neurons = np.zeros(MAX_DEPTH + 1, dtype=np.int)
#     unassign_all_nodes(G, labels)
#     label_graph(G, network, layer, neurons)

#     ##calculate error
#     mi_score = mutual_information(G, labels)
        
#     n, d = network.nodes(data=True)[0]
    
#     num_communities_detected = len(d['n'].nodes())
    
#     return mi_score, num_communities_detected

# def main_with_labels(params, gml_filename, labels, init=1, lam=10000):

#     return fitness(params['eta'], params['sigma'],
#                    params['e_sg'], params['e_en'], gml_filename, labels, init, lam)

# def main_no_labels(params, gml_filename, init=1, lam=10000):
    
#     G = nx.read_gml(gml_filename)
    
#     #start layer
#     layer = 0

#     #run ghsom algorithm
#     network, MQE = ghsom(G, lam, params['eta'],
#                          params['sigma'], np.inf, params['e_sg'], params['e_en'], init, layer)
    
        
#     n, d = network.nodes(data=True)[0]
    
#     return G, d['n']

def worker(q, networks):
    
    while True:
        
        #unpack first element of queue
        #conatains all the para,eters for GHSOM
        ID, G, lam, eta, sigma, e_0, e_sg, e_en, init, layer = q.get()
        
        #run GHSOM and return a network and MQE
        n, e = ghsom(ID, G, lam, eta, sigma, e_0, e_sg, e_en, init, layer, q)
        
        #append result to networks list
        networks.append((ID, n, e))
        
        #mark task as done
        q.task_done()

def main(params, gml_filename, init=1, lam=10000):
    
    #network
    G = nx.read_gml(gml_filename)
    
    ##list of networks
    networks = []
    
    #initilise worker queue
    q = Queue()
    
    #initialise threads
    for i in range(NUM_THREADS):
        
        t = Thread(target=worker, args=(q, networks))
        t.setDaemon(True)
        t.start()
        
    #add initial layer of ghsom to queue
    q.put(("00", G, lam, params["eta"], params["sigma"], np.inf, params["e_sg"], params["e_en"], init, 0))
    
    #finally wait until queue is empty and all tasks are done
    q.join()
    
    print "DONE"
    
    return G, networks


# In[2]:

params = {'eta': 0.0001,
         'sigma': 1,
          'e_sg': 0.3,
         'e_en': 0.3}


# In[ ]:

get_ipython().run_cell_magic(u'time', u'', u'G, networks = main(params=params, gml_filename="embedded_yeast_reactome.gml")')


# In[6]:

##TODO write networks to file 
for network in networks:
    
    ID, net, MQE = network
    
    print ID


# In[ ]:


