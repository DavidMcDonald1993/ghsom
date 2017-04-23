
# coding: utf-8

# In[1]:

import os
import pickle
import shutil
import Queue

import networkx as nx
import numpy as np
import pandas as pd

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

# for data in ["yeast_reactome", "yeast_uetz", "collins", "ccsb", "ito_core", "lc_multiple"]:
for data in ["yeast_reactome"]:
    
    print "data={}".format(data)
    print
    
    for e_sg in [0.3, 0.2]:
        
        print "e_sg={}".format(e_sg)
        print

        for e_en in [0.4, 0.3, 0.2, 0.1]:

            print "e_en={}".format(e_en)
            print

            os.chdir(root_dir)

            #ghsom parameters
            params = {'eta': 0.0001,
                     'sigma': 1,
                      'e_sg': e_sg,
                     'e_en': e_en}

            map_file = '{}_hierarchy_communities_{}_{}'.format(data, e_sg, e_en)

            if not os.path.isfile("{}.pkl".format(map_file)):

                #run ghsom and save output
                print "running GHSOM and saving to {}.pkl".format(map_file)
                G, map = ghsom_main(params, 'embedded_{}.gml'.format(data), init=1, lam=10000)
                print '\nnumber of communities detected: {}, saved map to {}'.format(len(map), map_file)
                save_obj((G, map), map_file)

            else:

                print "{}.pkl already exists, loading map".format(map_file)    
                #load output
                G, map = load_obj(map_file)

            #save results to file
            dir_name = "{}_hierarchy_communities_{}_{}".format(data, e_sg, e_en)
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
            print "wrote {}".format(all_genes_file)

            #map queue
            q = Queue.Queue()

            c = 1
            depth = 0
            q.put((c, depth, map))

            genes = G.nodes()
            gene_assignments = {k: v for k, v in zip(genes, 
                np.array([["" for j in range(10)] for i in range(len(genes))], dtype="S20"))}

            while not q.empty():

                map_id, depth, map = q.get()
                c = 1

                #shortest path matrix
                communities_in_this_map = np.array(["{}-{}".format(str(map_id).zfill(2), 
                                                                   str(i).zfill(2)) for i in range(c, c + len(map))])

                shortest_path = nx.floyd_warshall_numpy(map).astype(np.int)
    #             shortest_path = np.insert(shortest_path, 0, communities_in_this_map, axis=1)
                shortest_path_df = pd.DataFrame(shortest_path, index=communities_in_this_map)
                shortest_path_file = "{}_shortest_path.csv".format(str(map_id).zfill(2))
    #             np.savetxt(shortest_path_file, shortest_path, fmt='%i', delimiter=",")
                shortest_path_df.to_csv(shortest_path_file, index=True, header=False, sep=',')
                print 'wrote shortest path matrix and saved as {}'.format(shortest_path_file)



                #gene community assignments
                for n, d in map.nodes(data=True):

                    community = communities_in_this_map[c - 1]

                    for node in d['ls']:
                        gene_assignments[node][depth] = community

                    #add map to queue
                    m = d['n']

                    if not m == []:

                        q.put((community, depth + 1, m))   

                    c += 1

            #back to matrix
            assignment_matrix = np.array([v for k, v in gene_assignments.items()])
            #remove unnecessary columns
            mask = assignment_matrix != ""
            idx = mask.any(axis = 0)
            assignment_matrix = assignment_matrix[:,idx]
            assignment_matrix = np.insert(assignment_matrix, 0, "01", axis=1)

            assignment_matrix_file = "assignment_matrix.csv"
            np.savetxt(assignment_matrix_file, assignment_matrix, delimiter=",", fmt="%s")
            print "wrote assignment matrix and saved it as {}".format(assignment_matrix_file)
            print


# In[1]:

ID = "01-05"


# In[3]:

ID.count("-")


# In[ ]:



