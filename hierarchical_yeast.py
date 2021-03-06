
# coding: utf-8

# In[3]:

import os
import networkx as nx
import numpy as np
import pickle
import shutil
import Queue

from ghsom import main_no_labels as ghsom_main

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)
    
def add(d, key, value):
    if key in d:
        d[key].append(value)
    else:
        d.update({key : [value]})
    
root_dir = "/home/david/Documents/ghsom"

data = "yeast_uetz"
init = 1

for p in np.arange(0.1, 1, 0.1)[::-1]:
# for p in [0.05, 0.01]:
    
    print "p={}".format(p)
    
    os.chdir(root_dir)
    
    #ghsom parameters
    params = {'w': 0.0001,
             'eta': 0.0001,
             'sigma': 1,
              'e_sg': p,
             'e_en': p}
    
    map_file = '{}_hierarchy_communities_{}_{}'.format(data, p, init)
    
    if not os.path.isfile("{}.pkl".format(map_file)):
    
        #run ghsom and save output
        print "running GHSOM and saving to {}.pkl".format(map_file)
        G, map = ghsom_main(params, 'embedded_{}.gml'.format(data), init=init, lam=10000)
        print '\nnumber of communities detected: {}, saved map to {}'.format(len(map), map_file)
        save_obj((G, map), map_file)
    
    else:
        
        print "{}.pkl already exists, loading map".format(map_file)    
        #load output
        G, map = load_obj(map_file)

    #save results to file
    dir_name = "{}_hierarchy_communities_{}_{}".format(data, p, init)
    if os.path.isdir(dir_name):
        shutil.rmtree(dir_name)
        
    os.mkdir(dir_name)
    print 'made directory {}'.format(dir_name)

    os.chdir(dir_name)
    print "moved to {}".format(dir_name)
    
    #all genes
    all_genes_file = "all_genes.txt"
    with open(all_genes_file, 'w') as f:
        for n, d in G.nodes(data=True):
            f.write("{}\n".format(d['label']))
    print "written {}".format(all_genes_file)
    
    #map queue
    q = Queue.Queue()
    
    c = 1
    depth = 0
    q.put((c, depth, map))
    
    genes = G.nodes()
    gene_assignments = {k: v for k, v in zip(genes, -np.ones((len(genes), 10)))}
    

    while not q.empty():
        
        map_id, depth, map = q.get()
        
        #shortest path matrix
        communities_in_this_map = range(c + 1, c + 1 + len(map))
        shortest_path = nx.floyd_warshall_numpy(map).astype(np.int)
        shortest_path = np.insert(shortest_path, 0, communities_in_this_map, axis=1)
        shortest_path_file = "{}_shortest_path.csv".format(map_id)
        np.savetxt(shortest_path_file, shortest_path, fmt='%i', delimiter=",")
        print 'written shortest path matrix and saved as {}'.format(shortest_path_file)
        
        #gene community assignments
        for n, d in map.nodes(data=True):
            
            c += 1
            
            for node in d['ls']:
                gene_assignments[node][depth] = c
                
            #add map to queue
            m = d['n']
            
            if not m == []:
                
                q.put((c, depth + 1, m))
        
    #back to matrix
    assignment_matrix = np.array([v for k, v in gene_assignments.items()])
    #remove unnecessary columns
    mask = assignment_matrix != -1
    idx = mask.any(axis = 0)
    assignment_matrix = assignment_matrix[:,idx]
    assignment_matrix = np.insert(assignment_matrix, 0, 1, axis=1)
    
    assignment_matrix_file = "assignment_matrix.csv"
    np.savetxt(assignment_matrix_file, assignment_matrix, fmt='%i', delimiter=",")
    print "written assignment matrix and saved it as {}".format(assignment_matrix_file)
    print


# In[ ]:



