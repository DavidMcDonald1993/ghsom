
# coding: utf-8

# In[46]:

##imports
import numpy as np
import networkx as nx
import math
import matplotlib.pyplot as plt


# In[1]:

#function to generate real valued som for graph input
def initialise_network(X, num_neurons, w):
    
    #network will be a one dimensional list
    network = nx.Graph()
    
    #dimension of data in X
    d = len(X[0])
    
    for i in range(num_neurons):
        
        ##position
        network.add_node(i)
        
        ##weight    
        network.node[i]['v'] = 2 * (np.random.rand(d) - 0.5) * w

        ##list of closest nodes
        network.node[i]['ls'] = []

        ##error of neuron
        network.node[i]['e'] = 0
        
        ##som for neuron
        network.node[i]['n'] = []
        
        ## signal frequency
        network.node[i]['r'] = 0
        
        ##add edges
        for j in range(max(0, i-2), i):
            
            network.add_edge(i, j)
    
    #return network
    return network


# In[6]:

# function to train SOM on given graph
def train_network(X, network, num_epochs, eta_0, sigma_0, N):
    
    #initial learning rate
    eta = eta_0
    
    #initial neighbourhood size
    sigma = sigma_0
    
    #list if all patterns to visit
    training_patterns = [p for p in range(len(X))]
    
    for e in range(num_epochs):
        
        #shuffle nodes
        np.random.shuffle(training_patterns)
        
        # iterate through N nodes of graph
        for i in range(N):
            
            #data point to consider
            x = X[training_patterns[i]]
            
            #determine winning neuron
            win_neuron = winning_neuron(x, network)
            
            # update weights
            update_weights(x, network, win_neuron, eta, sigma)
            
        # drop neighbourhood
        sigma = sigma_0 * np.exp(-2 * sigma_0 * e / num_epochs);
        
#         ##decay signal frequency
#         for i in network.nodes():
#             network.node[i]['r'] -= 0.05 * network.node[i]['r'] 
    
#     return network


# In[43]:

# winning neuron
def winning_neuron(x, network):

    # minimum distance so far
    min_dist = math.inf
    win_neuron = []
    
    # iterate through network
    for i in network.nodes():
            
        #unpack network
        v = network.node[i]['v']
        
        if len(v) != len(x):
            print()
            print(v)
            print(x)
            print()

        #distance between input vector and neuron weight
        distance = np.linalg.norm(x - v)

        # if we have a new closest neuron
        if distance < min_dist:
            min_dist = distance
            win_neuron = i
            
    ##increment signal frequency
    network.node[win_neuron]['r'] += 1
    
    #return
    return win_neuron


# In[4]:

# function to update weights
def update_weights(x, network, win_neuron, eta, sigma):
    
    # iterate through all neurons in network
    for i in network.nodes():
        
        #unpack
        v = network.node[i]['v']

        #new v -- move along shortest path by move distance
        v += eta * neighbourhood(network, i, win_neuron, sigma) * (x - v)

        #save to network
        network.node[i]['v'] = v
        
        ##decay signal frequency
        network.node[i]['r'] -= 0.05 * network.node[i]['r'] 
    
#     return network


# In[3]:

# neighbourhood function
def neighbourhood(network, r, win_neuron, sigma):
    
    return np.exp(-(nx.shortest_path_length(network, r, win_neuron)) ** 2 / (2 * sigma ** 2))


# In[1]:

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
        min_distance = math.inf
        
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
    
    #return assignment list
#     return network
            


# In[15]:

##function to return lattice grid of errors
def update_errors(network):
    
    #mean network error
    mqe = 0;
    
    #number of nodes
    num_neurons = nx.number_of_nodes(network)
    
    #iterate over all neurons and average distance
    for i in network.nodes():
            
        #unpack network
        e = network.node[i]['e']
        ls = network.node[i]['ls']

        if len(ls) == 0:
            continue

        #divide by len(ls) for mean
        e /= len(ls)

        #sum total errors
        mqe += e
        
        #save error to network
        network.node[i]['e'] = e      
    
    #mean
    mqe /= num_neurons
    
    #return network
#     return network, mqe
    return mqe


# In[26]:

##function to identify neuron with greatest error
def identify_error_unit(network):
    
    #first node
    n = network.nodes()[0]
    
    #initial value for maximum error found
    max_e = network.node[n]['e']
    
    #initial index to return
    error_node = n
    
    for i in network.nodes():
        
        #unpack
        e = network.node[i]['e']
        
        #check if this unit has greater error than maximum 
        if e > max_e:
            
            max_e = e
            error_node = i
            
    #return id of unit with maximum error
    return error_node


# In[20]:

##function to expand som using given error unit
def expand_network(network, error_unit):
    
    #identify neighbour pointing furthest away
    error_unit_neighbours = network.neighbors(error_unit)
    
    #id of new node
    id = nx.number_of_nodes(network)
    
    #v of error unit
    ve = network.node[error_unit]['v']
    
    #dimension
    d = len(ve)
    
    ##
    if len(error_unit_neighbours) == 0:
        ##random position
        
        ##position
        network.add_node(id)
        
        ##weight    
        v = 2 * (np.random.rand(d) - 0.5) * 1e-2
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
        
        ## signal frequency
        r = 0
        network.node[id]['r'] = r
        
        ##add edge
        network.add_edge(error_unit, id)
        
    elif len(error_unit_neighbours) == 1:
        
        #neighbour
        neighbour = error_unit_neighbours[0]
        
        #v of neighbour
        vn = network.node[neighbour]['v']
        
        ##position
        network.add_node(id)
        
        ##weight    
        v = (ve + vn) / 2
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
        
        ## signal frequency
        r = 0
        network.node[id]['r'] = r
        
        ##add edge
        network.add_edge(error_unit, id)
        network.add_edge(neighbour, id)
        
    else:

        #neighbour id
        n1 = furthest_neuron(network, error_unit, error_unit_neighbours)
        
        ##v of n1
        v_n1 = network.node[n1]['v']
            
        #now we have identified neighbour pointing furthest away in input space
        #take mean and produce new neuron
        
        ##must find mutual neighbours
        neighbour_neighbours = network.neighbors(n1)
        
        ##mutual neighbours
        mutual_neighbours = [n for n in error_unit_neighbours if n in neighbour_neighbours]

        ##second furthest node
        n2 = furthest_neuron(network, error_unit, mutual_neighbours)
        
        #v of n2
        v_n2 = network.node[n2]['v']
        
        ##position
        network.add_node(id)
        
        ##weight    
        v = (ve + v_n1 + v_n2) / 3 
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
        
        ## signal frequency
        r = 0
        network.node[id]['r'] = r
        
        #remove edge from n1 and n2
        network.remove_edge(n1, n2)
        
        #connect new node to all nodes that are connected to both neighbours (including error unit)
        for neuron in network.nodes():
            if network.has_edge(n1, neuron) and network.has_edge(n2, neuron):
                ##add edges
                network.add_edge(neuron, id)
        
        network.add_edge(n1, id)
        network.add_edge(n2, id)
    
    #return network
#     return network


# In[42]:

##function to expand som using given error unit
def expand_network2(G, network, error_unit, layer):
    
    #id of new node
    id = nx.number_of_nodes(network)
    
    #v of error unit
    ve = network.node[error_unit]['v']
    
    #dimension
    d = len(ve)
    
    ##position
    network.add_node(id)

    ##weight    
    network.node[id]['v'] = ve

    ##list of closest nodes
    network.node[id]['ls'] = []

    ##error of neuron
    network.node[id]['e'] = 0

    ##som for neuron
    network.node[id]['n'] = []

    ## signal frequency
    network.node[id]['r'] = 0
    
    #point error unit to furthest datapoint
    network.node[error_unit]['v'] = furthest_datapoint(G, network, error_unit, layer)
    
    #identify neighbour pointing furthest away
    error_unit_neighbours = network.neighbors(error_unit)
    
    print('furthest datapoint',network.node[error_unit]['v'])
    
    ##
    if len(error_unit_neighbours) == 0:
        
        ##add edge
        network.add_edge(error_unit, id)
        
    elif len(error_unit_neighbours) == 1:
        
        #neighbour
        neighbour = error_unit_neighbours[0]
        
        ##add edge
        network.add_edge(error_unit, id)
        network.add_edge(neighbour, id)
        
    else:

        #neighbour id
        n1 = furthest_neuron(network, error_unit, error_unit_neighbours)
            
        #now we have identified neighbour pointing furthest away in input space
        #take mean and produce new neuron
        
        ##must find mutual neighbours
        neighbour_neighbours = network.neighbors(n1)
        
        ##mutual neighbours
        mutual_neighbours = [n for n in error_unit_neighbours if n in neighbour_neighbours]

        ##second furthest node
        n2 = furthest_neuron(network, error_unit, mutual_neighbours)
        
        #remove edge from n1 and n2
        network.remove_edge(n1, n2)
        
        #connect new node to all nodes that are connected to both neighbours (including error unit)
        for neuron in network.nodes():
            if network.has_edge(n1, neuron) and network.has_edge(n2, neuron):
                ##add edges
                network.add_edge(neuron, id)
        
        network.add_edge(n1, id)
        network.add_edge(n2, id)
    
    #return network
#     return network


# In[41]:

##find vector of furthest datapoint
def furthest_datapoint(G, network, n, layer):
    
    v = network.node[n]['v']
    ls = network.node[n]['ls']
    
    if len(ls) == 0:
        print('empty list',network.node[n]['e'])
    
    max_dist = -math.inf
    furthest_x = []
    
    for node in ls:
        
#         x = G.node[node]['embedding'+str(layer)]
        x = G.node[node]['embedding0']
    
        d = np.linalg.norm(x - v)
        
        if d > max_dist:
            
            max_dist = d
            furthest_x = x
            
    return furthest_x


# In[42]:

##function to find neuron pointing furthest away in list
def furthest_neuron(network, error_unit, ls):
    
    vi = network.node[error_unit]['v']
    max_dist = -math.inf
    
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


# In[22]:

##delete a neuron 
def delete_neurons(network):
    
    deleted = False
    
#     print('number of neurons:',len(network))
    
    for n in network.nodes():
        
        r = network.node[n]['r']
        
#         print('r',r)
        
        if r < 0.09:
            
            print('deleting neuron',n)

            neighbours = network.neighbors(n)

            num_neighbours = len(neighbours)

            for n1 in range(num_neighbours):

                for n2 in range(n1+1,num_neighbours):

                    network.add_edge(neighbours[n1], neighbours[n2])

            network.remove_node(n)
            
            deleted = True
    
    return deleted


# In[5]:

##function to visualise graph
def visualise_graph(G, colours, layer):
        
    ## create new figure for graph plot
    fig, ax = plt.subplots()
    
    # graph layout
    pos = nx.spring_layout(G)
    
    #attributes in this graph
    attributes = np.unique([v for k,v in nx.get_node_attributes(G, 'community'+str(layer)).items()])

    # draw nodes -- colouring by cluster
    for i in range(min(len(colours), len(attributes))):
       
        node_list = [n for n in G.nodes() if G.node[n]['community'+str(layer)] == attributes[i]]
        colour = [colours[i] for n in range(len(node_list))]
        
        nx.draw_networkx_nodes(G, pos, nodelist=node_list, node_color=colour)
        
    #draw edges
    nx.draw_networkx_edges(G, pos)

    # draw labels
    nx.draw_networkx_labels(G, pos, )
    
    #title of plot
    plt.title('Nodes coloured by cluster, layer: '+str(layer))

    #show plot
    plt.show()


# In[4]:

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

