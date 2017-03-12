
# coding: utf-8

# In[2]:

# import karate graph
import networkx as nx

# import karate club graph
G=nx.karate_club_graph()

#distance matrix
dist = nx.all_pairs_shortest_path_length(G)

# number of nodes in G
num_nodes = len(G.nodes())

# lattice size
lattice_size_x = 2
lattice_size_y = 1

#initial neighborhood size
sigma_0 = 3

#initial learning rate
eta_0 = 0.01

#number of training epochs
num_epochs = 1000

# num_inputs to train on per epoch
N = num_nodes


# In[3]:

import som_functions as som

## MAIN PROGRAM

network = som.initialise_network(G, lattice_size_x, lattice_size_y)

print("initialised network")

network = som.train_network(G, network, dist, num_epochs, eta_0, sigma_0, N)

print("trained network")


# In[4]:

import som_functions as som


#visualise based on som clusters
som.visualise_graph(G, network, dist)


# In[3]:

## visualise karate graph actual classes

#draw with networkx draw
import matplotlib.pyplot as plt

# get class labels
clubs = nx.get_node_attributes(G, 'club')

# divide nodes into clubs
mr_hi = [k for k,v in clubs.items() if v == 'Mr. Hi']
officer = [k for k,v in clubs.items() if v == 'Officer']

# graph layout
pos = nx.fruchterman_reingold_layout(G)

# draw nodes -- colouring by club
nx.draw_networkx_nodes(G, pos, nodelist = mr_hi, node_color = 'r')
nx.draw_networkx_nodes(G, pos, nodelist = officer, node_color = 'b')

#draw edges
nx.draw_networkx_edges(G, pos)

# draw labels
nx.draw_networkx_labels(G, pos)

#show plot
plt.show()

