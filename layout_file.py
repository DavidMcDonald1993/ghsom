
# coding: utf-8

# In[19]:

def write_layout_file(G, label, filename):
    
    with open('{}.layout'.format(filename), 'w') as f:
        
        for e in G.edges_iter():
            f.write('\"{}\"\t\"{}\"\n'.format(e[0], e[1])) 

        for n in G.nodes():
            if label in G.node[n]:
                f.write('//NODECLASS\t\"{}\"\t\"{}\"\n'.format(n, G.node[n][label]))


# In[20]:

import networkx as nx

G = nx.read_gml('dolphins_labelled.gml')

write_layout_file(G, 'group', 'dolphin')


# In[2]:

import networkx as nx

G = nx.read_gml("embedded_yeast_union.gml")

len(G)


# In[7]:

import os
import networkx as nx
import numpy as np
from spearmint_ghsom import main_no_labels as ghsom_main
import pickle
import shutil

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)
    
map_file = '{}_communities_{}'.format(data, p, init)

G, map = load_obj("yeast_union_communities_0.6_1")


# In[2]:

data = "union"
p = 0.6
init = 1

dir_name = "{}_communities_{}_{}".format(data, p, init)

os.chdir(dir_name)


# In[4]:




# In[ ]:



