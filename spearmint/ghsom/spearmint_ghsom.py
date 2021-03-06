
import numpy as np
import sys
import networkx as nx
import numpy as np
import networkx as nx
import sklearn.metrics as met
import sklearn.manifold as man

MAX_DEPTH = 2

#function to generate real valued som for graph input
def initialise_network(X, num_neurons, w):
    
    #network will be a one dimensional list
    network = nx.Graph()
    
    #number of data points
    num_nodes = len(X) 
    
    #dimension of data in X
    d = len(X[0])
    
    for i in range(num_neurons):
        
        ##position
        network.add_node(i)
        
        ##weight    
#         network.node[i]['v'] = 2 * (np.random.rand(d) - 0.5) * w
        r = np.random.randint(num_nodes)
        network.node[i]['v'] = X[r]

        ##list of closest nodes
        network.node[i]['ls'] = []

        ##error of neuron
        network.node[i]['e'] = 0
        
        ##som for neuron
        network.node[i]['n'] = []
        
        # signal frequency
        network.node[i]['r'] = 0
        
        ##add edges
        for j in range(max(0, i-2), i):
            
            network.add_edge(i, j)
    
    #return network
    return network


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

# winning neuron
def winning_neuron(x, network):

    # minimum distance so far
    min_dist = float("inf")
    win_neuron = []
    
    # iterate through network
    for i in network.nodes():
            
        #unpack network
        v = network.node[i]['v']

        #distance between input vector and neuron weight
        distance = np.linalg.norm(x - v)

        # if we have a new closest neuron
        if distance < min_dist:
            min_dist = distance
            win_neuron = i
            
    #increment signal frequency
#     network.node[win_neuron]['r'] += 1
    
    #return
    return win_neuron

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

# neighbourhood function
def neighbourhood(network, r, win_neuron, sigma):
    
    return np.exp(-(nx.shortest_path_length(network, r, win_neuron)) ** 2 / (2 * sigma ** 2))

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
        min_distance = float("inf")
        
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
def update_errors(network):
    
    #mean network error
    mqe = 0;
    
    nodes_to_remove = []
    
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
            
    network.remove_nodes_from(nodes_to_remove)        
    
    #mean
    mqe /= nx.number_of_nodes(network)
    
    return mqe

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
        
#         ## signal frequency
#         r = 0
#         network.node[id]['r'] = r
        
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
        
#         ## signal frequency
#         r = 0
#         network.node[id]['r'] = r
        
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
        
#         ## signal frequency
#         r = 0
#         network.node[id]['r'] = r
        
        #remove edge from n1 and n2
        network.remove_edge(n1, n2)
        
        #connect new node to all nodes that are connected to both neighbours (including error unit)
        for neuron in network.nodes():
            if network.has_edge(n1, neuron) and network.has_edge(n2, neuron):
                ##add edges
                network.add_edge(neuron, id)
        
        network.add_edge(n1, id)
        network.add_edge(n2, id)

##function to find neuron pointing furthest away in list
def furthest_neuron(network, error_unit, ls):
    
    vi = network.node[error_unit]['v']
    max_dist = -float("inf")
    
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

##GHSOM algorithm
def ghsom(G, lam, w, eta, sigma, e_0, e_sg, e_en, layer):
    
    #embedding
    X = get_embedding(G)
    
    #number of nodes in G
    num_nodes = nx.number_of_nodes(G)
    
    ##number of training patterns to visit
    N = min(num_nodes, 100)
    
    #create som for this neuron
    network = initialise_network(X, 1, w)
    
    ##inital training phase
    
    #train for lam epochs
    train_network(X, network, lam, eta, sigma, N)

    #classify nodes
    assign_nodes(G, X, network, layer)

    #calculate mean network error
    MQE = update_errors(network)
    
    ##som growth phase
    #repeat until error is low enough
    while MQE > e_sg * e_0:
    
        #find neuron with greatest error
        error_unit = identify_error_unit(network)
        
        #expand network
        expand_network(network, error_unit)
        
        #train for l epochs
        train_network(X, network, lam, eta, sigma, N)

        #classify nodes
        assign_nodes(G, X, network, layer)

        #calculate mean network error
        MQE = update_errors(network)
    
    #recalculate error after neuron expansion
    MQE = 0
    
    ##neuron expansion phase
    #iterate thorugh all neruons and find neurons with error great enough to expand
    for i in range(len(network)):
        
        #unpack
        ls = network.node[network.nodes()[i]]['ls']
        e = network.node[network.nodes()[i]]['e']
        
        #check error
        if (e > e_en * e_0 and layer < MAX_DEPTH) or e_0 == float("inf"):

            if e_0 == float("inf"):
                e_0 = e
        
            #subgraph
            H = G.subgraph(ls)
            
            #recursively run algorithm to create new network for subgraph of this neurons nodes
            n, e = ghsom(H, lam, w, eta, sigma, e_0, e_sg, e_en, layer + 1)
            
            #repack
            network.node[network.nodes()[i]]['e'] = e
            network.node[network.nodes()[i]]['n'] = n
            
        #increase overall network error
        MQE += e
    
    #mean MQE
    MQE /= nx.number_of_nodes(network)
    
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

def modularity(network):
    
    #modularity
    Q = 0

#     for neuron in network:
        
    
    return Q

###evaluate fitness
def fitness(w, eta, sigma, e_sg, e_en, gml_filename, labels, lam):

    G = nx.read_gml(gml_filename)
    labels = labels.split(',')

    #start layer
    layer = 0

    #run ghsom algorithm
    network, MQE = ghsom(G, lam, w, eta, sigma, float("inf"), e_sg, e_en, layer)

    #label graph
    neurons = np.zeros(MAX_DEPTH + 1, dtype=np.int)
    unassign_all_nodes(G, labels)
    label_graph(G, network, layer, neurons)

    ##calculate error
    mi_score = mutual_information(G, labels)

    #add to list of scores
    mi_scores = mi_score
        
    n, d = network.nodes(data=True)[0]
    
    num_communities_detected = len(d['n'].nodes())
    
    return mi_scores, num_communities_detected

def main(job_id, params, gml_filename, labels, lam):

    print 'running job id: {}'.format(job_id)

    return 1 - fitness(0.0001, 0.0001, 1,
                   params['e_sg'], 0.8, gml_filename, labels, lam)[0]

def main_no_labels(params, gml_filename, lam):
    
    G = nx.read_gml(gml_filename)
    
    #start layer
    layer = 0

    #run ghsom algorithm
    network, MQE = ghsom(G, lam, params['w'], params['eta'],
                         params['sigma'], float("inf"), params['e_sg'], params['e_en'], layer)
        
    n, d = network.nodes(data=True)[0]
    
    return G, d['n']

