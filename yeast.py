
# coding: utf-8

# In[24]:

import numpy as np

filename = "HI-II-14.tsv"

lines = []

for line in open(filename):
    lines.append(line.rstrip().split("\t"))

with open("HI-II-14.txt", "w") as f:
    first = True
    for line in lines:
        if first:
            first=False
            continue
        f.write("{} {}\n".format(line[0], line[2]))


# In[27]:

import networkx as nx
G = nx.read_edgelist("HI-II-14.txt")

print len(G.nodes())


# In[1]:

from save_embedded_graph27 import main as embed_main

embed_main('HI-II-14.txt', 'embedded_hi_ii_14.gml')


# In[10]:

from spearmint_ghsom import main_no_labels as ghsom_main
import pickle
import os

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)
    
os.chdir("/home/david/Documents/ghsom")

p = 0.8

#ghsom parameters
params = {'w': 0.0001,
         'eta': 0.0001,
         'sigma': 1,
          'e_sg': p,
         'e_en': 0.8}

G, map = ghsom_main(params, 'embedded_hi_ii_14.gml')

print 'number of communities detected: {}'.format(len(map))
save_obj((G, map), 'HI_II_communities_{}'.format(p))

print 'done'


# In[11]:

import os

os.chdir("/home/david/Documents/ghsom")

G, map = load_obj('HI_II_communities_{}'.format(p))
print 'num communities: {}'.format(len(map))


# In[12]:

min_nodes = 10

##remove neurons with no assigned nodes
for n, d in map.nodes(data=True):
    if len(d['ls']) < min_nodes:
        map.remove_node(n)
        print 'removed node {}'.format(n)


# In[13]:

import os
import networkx as nx
import numpy as np
##save to communities directory
os.chdir("/home/david/Documents/ghsom")

dir_name = "hi_communities_08"

if not os.path.isdir(dir_name):
    os.mkdir(dir_name)
    print 'made directory {}'.format(dir_name)
    
os.chdir(dir_name)

shortest_path = nx.floyd_warshall_numpy(map).astype(np.int)
np.savetxt("shortest_path.csv", shortest_path, fmt='%i', delimiter=",")
print 'written shortest path matrix'

c = 0
for n, d in map.nodes(data=True):
    ls = d['ls']
    with open('community_{}.txt'.format(c),'w') as f:
        for l in ls:
            f.write('{}\n'.format(l))
    print 'written community_{}.txt'.format(c)
    c += 1


# In[ ]:



