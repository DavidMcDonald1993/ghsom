
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


# In[4]:

import numpy as np

nmi_scores = np.genfromtxt('nmi_scores.csv', delimiter=',')

print nmi_scores

print np.mean(nmi_scores, axis=0)


# In[ ]:



