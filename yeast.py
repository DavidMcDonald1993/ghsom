
# coding: utf-8

# In[15]:

import numpy as np

os.chdir("/home/david/Documents/ghsom")

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


# In[1]:

import networkx as nx
G = nx.read_edgelist("Uetz_screen.txt")

print len(G.nodes())


# In[1]:

from save_embedded_graph27 import main as embed_main

embed_main('HI-II-14.txt', 'embedded_hi_ii_14.gml')


# In[ ]:

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

init = 1
p = 0.3

#ghsom parameters
params = {'w': 0.0001,
         'eta': 0.001,
         'sigma': 1,
          'e_sg': p,
         'e_en': 10}

# G, map = ghsom_main(params, 'embedded_hi_ii_14.gml')
G, map = ghsom_main(params, 'embedded_yeast_union.gml', init=init, lam=1000)
# G, map = ghsom_main(params, 'embedded_yeast_uetz.gml', init=init, lam=1000)

print '\nnumber of communities detected: {}'.format(len(map))
# save_obj((G, map), 'HI_II_communities_{}'.format(p))
# save_obj((G, map), 'yeast_union_communities_{}'.format(p))
save_obj((G, map), 'yeast_uetz_communities_{}_{}'.format(p, init))

print 'done'


# In[2]:

import os

os.chdir("/home/david/Documents/ghsom")

# G, map = load_obj('HI_II_communities_{}'.format(p))
# G, map = load_obj('yeast_union_communities_{}'.format(p))
G, map = load_obj('yeast_uetz_communities_{}_{}'.format(p, init))
print 'number of communities detected: {}'.format(len(map))


# In[12]:

min_nodes = 10

##remove neurons with no assigned nodes
for n, d in map.nodes(data=True):
    if len(d['ls']) < min_nodes:
        map.remove_node(n)
        print 'removed node {}'.format(n)


# In[5]:

import os
import networkx as nx
import numpy as np
##save to communities directory
os.chdir("/home/david/Documents/ghsom")

# dir_name = "union_communities_08"
dir_name = "uetz_communities_{}_{}".format(p, init)

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
    
root_dir = "/home/david/Documents/ghsom"

data = "hi_ii_14"
init = 1

for p in np.arange(0.1, 1, 0.1)[::-1]:
    
    print "p={}".format(p)
    
    os.chdir(root_dir)
    
    #ghsom parameters
    params = {'w': 0.0001,
             'eta': 0.001,
             'sigma': 1,
              'e_sg': p,
             'e_en': 10}
    
    map_file = '{}_communities_{}'.format(data, p, init)
    
    if not os.path.isfile("{}.pkl".format(map_file)):
    
        #run ghsom and save output
        print "running GHSOM and saving to {}.pkl".format(map_file)
        G, map = ghsom_main(params, 'embedded_{}.gml'.format(data), init=init, lam=1000)
        print '\nnumber of communities detected: {}, saved map to {}'.format(len(map), map_file)
        save_obj((G, map), map_file)
    
    else:
        
        print "{}.pkl already exists, loading map".format(map_file)    
        #load output
        G, map = load_obj(map_file)

    #save results to file
    dir_name = "{}_communities_{}_{}".format(data, p, init)
    if os.path.isdir(dir_name):
        shutil.rmtree(dir_name)
        print "deleted directory {}".format(dir_name)
    
    os.mkdir(dir_name)
    print 'made directory {}'.format(dir_name)

    os.chdir(dir_name)
    print "moved to {}".format(dir_name)
    
    #all genes
    all_genes_file = "all_genes.txt"
    with open(all_genes_file, 'w') as f:
        for n in G.nodes():
            f.write("{}\n".format(n))
    print "written {}".format(all_genes_file)
    
    #save shortest path matrix
    shortest_path = nx.floyd_warshall_numpy(map).astype(np.int)
    np.savetxt("shortest_path.csv", shortest_path, fmt='%i', delimiter=",")
    print 'written shortest path matrix'
    
    #save communities to file
    c = 0
    for n, d in map.nodes(data=True):
        ls = d['ls']
        with open('community_{}.txt'.format(c),'w') as f:
            for l in ls:
                f.write('{}\n'.format(l))
        print 'written community_{}.txt'.format(c)
        c += 1
    print


# In[ ]:



