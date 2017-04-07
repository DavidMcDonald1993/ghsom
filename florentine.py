
# coding: utf-8

# In[9]:

import numpy as np
d = 10000
a = np.random.rand(d, d)


# In[10]:

get_ipython().run_cell_magic(u'time', u'', u'b = np.zeros((d, d))\n\nfor i in range(d):\n    for j in range(d):\n        b[i,j] = a[i,j] **2')


# In[ ]:

b[:,:10]


# In[ ]:

get_ipython().run_cell_magic(u'time', u'', u'c = np.array([[j ** 2 for j in i] for i in a])')


# In[ ]:

c[:,:10]


# In[1]:

import networkx as nx

G = nx.florentine_families_graph()


# In[3]:

G.nodes(data=True)


# In[4]:

G.edges()


# In[8]:

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

data = "florentine_families"
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
    all_genes_file = "all_families.txt"
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


# In[1]:

import networkx as nx

G = nx.read_gml("embedded_florentine_families.gml")


# In[5]:

from spearmint_ghsom import get_embedding

X = get_embedding(G)


# In[6]:

X


# In[ ]:



